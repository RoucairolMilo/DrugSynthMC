use std::fs::{self, File};
use crate::models::SMILESgen::{State, Move};
use crate::tools::calc::softmaxChoice;
use crate::tools::resultSaver::writeLine;
use std::time::Instant;

//c'est en classe uniquement pour avoir accès à la variable globale
pub struct NMCS{
    pub best_yet : f64,
    pub timeout : f64,
    pub registerName : String,
    pub start_time : Instant,
}

impl NMCS{
    pub fn new() -> Self {
        Self{ start_time : Instant::now(),best_yet: 0.0, timeout : -1.0, registerName : String::new() }
    }

    pub fn playout(&mut self, mut st: State, heuristic_w : f64) -> State{
        let mut best_state: State = st.clone(); //State::new();

        let mut best_state_score = 0.0;
        if State::CONSIDER_NON_TERM || st.terminal() {
            best_state_score = best_state.score();
        }

        while !st.terminal() {
            let moves = st.legal_moves();
            if moves.len() == 0 { break }

            let mut i = ((moves.len() as f64)*rand::random::<f64>()) as usize;
            if heuristic_w != 0.0 {
                let mut weights = Vec::new();
                for &m in &moves{
                    weights.push(heuristic_w*st.heuristic(m));
                }
                i = softmaxChoice(weights);
            }
            let mv = moves[i];


            /*
            //greedy myope
            let mut mv = moves[0];
            let mut best = 0.0;
            for m in moves {
                let heuri = st.heuristic(m);
                if heuri > best {
                    best = heuri;
                    mv = m;
                }
            }

             */



            st.play(mv);

            if State::CONSIDER_NON_TERM {
                let sc = st.score();
                if sc > best_state_score {
                    best_state_score = sc;
                    best_state = st.clone();
                }
            }
        }

        //println!("playout done");
        if State::CONSIDER_NON_TERM{
            return best_state;
        }
        return st;
    }

    pub fn nmcs(&mut self, mut st: State, n : i8, heuristic_w : f64, verbose : bool) -> State{

        let mut best_state: State = st.clone(); //State::new();
        let mut best_state_score = -1.0; //best_state.score();

        while !st.terminal(){
            let moves = st.legal_moves();
            if moves.len() == 0 { break }
            for &mv in &moves{
                if self.start_time.elapsed().as_secs_f64() > self.timeout && self.timeout > 0.0 {return best_state;}

                let mut new_st = st.clone();
                new_st.play(mv);
                if n <= 1 {
                    new_st = self.playout(new_st, heuristic_w);
                }else{
                    new_st = self.nmcs(new_st, n-1, heuristic_w, verbose);
                }
                let new_st_score = new_st.score();

                //println!("le score : {}", new_st_score);


                //HP model et SMILEGEN

                if new_st.reached_best_score {
                    //on arrête tout
                    if verbose{ println!("reached best score !!!");}
                    return new_st
                }


                if new_st_score > best_state_score {
                    best_state = new_st;
                    best_state_score = new_st_score;
                    if best_state_score > self.best_yet {
                        self.best_yet = best_state_score;
                        let elapsed = self.start_time.elapsed().as_secs_f64();
                        writeLine(elapsed.to_string() + " " + &*best_state_score.to_string()+ "\n", self.registerName.clone());
                        if verbose {
                            println!("NMCS best score yet : {} after {}", best_state_score, elapsed);
                            //println!("{}", best_state.adj_mat);
                        }
                    }
                }
            }
            if State::CONSIDER_NON_TERM {
                if best_state.seq.len() == st.seq.len() {
                    break
                }
            }



            st.play(best_state.seq[st.seq.len()]);
        }


        if State::CONSIDER_NON_TERM{
            return best_state;
        }
        return st;
    }
}

pub fn launch_nmcs(init_st : State, level : i8, heuristic_w : f64, verbose : bool, timeout : f64, registerName : String) -> State{
    let mut expe = NMCS::new();
    expe.timeout = timeout;
    expe.registerName = registerName;

    //let mut st = State::new();
    let mut st = expe.nmcs(init_st, level, heuristic_w, verbose);
    return st;
}
