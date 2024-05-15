use std::thread;
use std::thread::JoinHandle;
use std::time::Instant;
use crate::models::SMILESgen::{State, Move};
use crate::tools::calc::softmaxChoice;
use crate::tools::resultSaver::writeLine;

//on utilise la moyenne au lieu du max
pub struct lazyNMCSv2{
    pub start_time : Instant,
    pub best_yet : f64,
    pub tresh : Vec<mean>,
    pub timeout : f64,
    pub sliding : i32,
    pub registerName : String,

    pub pruned : i32,
    pub spared : i32,
}

pub struct mean{
    pub val : f64,
    pub n : f64
}

pub fn playout(mut st: State, heuristic_w : f64) -> State{

    let mut best_state: State = st.clone(); //State::new();
    let mut best_state_score = best_state.score();

    while !st.terminal() {
        let moves = st.legal_moves();
        if moves.len() == 0 {return st;}
        let mut i = ((moves.len() as f64)*rand::random::<f64>()) as usize;
        if heuristic_w != 0.0 {
            let mut weights = Vec::new();
            for &m in &moves{
                weights.push(heuristic_w*st.heuristic(m));
            }
            i = softmaxChoice(weights);
        }
        let mv = moves[i];
        st.play(mv);

        if State::CONSIDER_NON_TERM {
            let sc = st.score();
            if sc > best_state_score {
                best_state_score = sc;
                best_state = st.clone();
            }
        }
    }

    if State::CONSIDER_NON_TERM{
        return best_state;
    }
    return st;
}

impl lazyNMCSv2{
    pub fn new() -> Self {
        Self{start_time : Instant::now(),
            best_yet: 0.0, tresh : vec![], timeout : 150.0, sliding : 0, registerName : String::new(), pruned : 0, spared : 0 }
    }

    pub fn lazy_nmcs_v2(&mut self, mut st: State, n : i8, p : i8, ratio : f64, heuristic_w : f64, verbose : bool) -> State{
        let mut best_state: State = st.clone(); //State::new();
        let mut best_state_score = -1.0; //best_state.score();
        while !st.terminal(){
            let moves = st.legal_moves();
            println!("{}", self.start_time.elapsed().as_secs_f64());
            if moves.len() == 0 || (self.start_time.elapsed().as_secs_f64() > self.timeout && self.timeout > 0.0) {return best_state;}

            for &mv in &moves{
                if self.start_time.elapsed().as_secs_f64() > self.timeout && self.timeout > 0.0 {return best_state;}

                let mut new_st = st.clone();
                new_st.play(mv);


                if n <= 1 {
                    new_st = playout(new_st, heuristic_w);
                }else{
                    let mut handles : Vec<JoinHandle<f64>> = Vec::new();
                    let mut eval_Scores = 0.0;
                    for _ in 0..p{
                        let mut eval_st = new_st.clone();
                        let handle = thread::spawn( move || {
                            eval_st = playout(eval_st, heuristic_w);
                            eval_st.score()
                        });
                        handles.push(handle);
                    }

                    for h in handles {
                        let res = h.join().unwrap();
                        eval_Scores += res;
                    }

                    let evalScore :f64 = eval_Scores/(p as f64);

                    //evalS.push(evalScore);

                    if self.tresh.len() < st.seq.len()+1{
                        self.tresh.push(mean{ val: 0.0, n: 0.0 });
                    }

                    if self.sliding == 0{
                        self.tresh[st.seq.len()].val = (self.tresh[st.seq.len()].val *  self.tresh[st.seq.len()].n + evalScore)/(self.tresh[st.seq.len()].n +1.0);
                    }else{
                        let mut s = self.sliding as f64;
                        if s > self.tresh[st.seq.len()].n{
                            s = self.tresh[st.seq.len()].n;
                        }
                        self.tresh[st.seq.len()].val = self.tresh[st.seq.len()].val*s /(s+1.0) + evalScore/(s+1.0);
                    }
                    self.tresh[st.seq.len()].n += 1.0;

                    //println!("{}", self.tresh[st.seq.len()].val);
                    if evalScore < ratio*self.tresh[st.seq.len()].val {
                        //println!("prunning!");

                        self.pruned += 1;
                        new_st = playout(new_st.clone(), heuristic_w);
                    }else{
                        //println!("spared!");
                        self.spared += 1;
                        new_st = self.lazy_nmcs_v2(new_st, n-1, p, ratio, heuristic_w, verbose);
                    }

                    //new_st = self.lazy_nmcs(new_st, n-1, p, heuristic_w, verbose);
                    //playS.push(new_st.score())

                    //println!("monomere {}, score Playout {}, score NMCS {} {}", st.molecule.len(), new_st0.score(), n-1, new_st.score());
                }
                let new_st_score = new_st.score();

                /*
                if new_st.reached_best_score {
                    //on arrête tout
                    if verbose{ println!("reached best score !!!");}
                    return new_st
                }

                 */


                if new_st_score > best_state_score {
                    best_state = new_st;
                    best_state_score = new_st_score;
                    if best_state_score > self.best_yet {
                        println!(" pruned : {}  spared : {}", self.pruned, self.spared);
                        self.best_yet = best_state_score;
                        let elapsed = self.start_time.elapsed().as_secs_f64();
                        writeLine(elapsed.to_string() + " " + &*best_state_score.to_string()+ "\n", self.registerName.clone());
                        if verbose {println!("best score yet : {} after {}", best_state_score, elapsed);}


                    }
                }
            }

            /*
            for i in 0..evalS.len(){
                print!("{} ", evalS[i]/playS[i]);

            }
            if n > 1 { println!(" ");}


             */
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

pub fn launch_lazy_nmcs_v2(init_stat : State, level : i8, p : i8, ratio : f64, heuristic_w : f64, timeout : f64, sliding : i32, verbose : bool, registerName : String) -> State{
    let mut expe = lazyNMCSv2::new();
    expe.timeout = timeout;
    expe.sliding = sliding;
    expe.registerName = registerName;
    //let mut st = State::new();

    let mut st = expe.lazy_nmcs_v2(init_stat, level, p, ratio, heuristic_w, verbose);

    return st;
}
