use std::thread;
use std::thread::JoinHandle;
use std::time::Instant;
//use crate::models::conjectures::GenerateGraph::{State};
//use crate::models::HPmodel::{State};
use crate::models::Network::{State};
use crate::tools::calc::softmaxChoice;
use crate::tools::resultSaver::writeLine;

//on utilise la moyenne au lieu du max
pub struct lazyNMCSv3{
    pub start_time : Instant,
    pub best_yet : f64,
    pub tresh : Vec<mean>,
    pub treshmax : Vec<f64>,
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

pub struct WS{
    pub w :f64,
    pub s : State
}

pub fn playout(mut st: State, heuristic_w : f64) -> State{

    let mut best_state: State = st.clone(); //State::new();
    let mut best_state_score = best_state.score();

    while !st.terminal() {

        let moves = st.legal_moves();

        if moves.len() == 0 {return st;}
        //base
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
            //let mesure = Instant::now();
            let sc = st.score();
            //println!("{}", mesure.elapsed().as_secs_f64());

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

impl lazyNMCSv3{
    pub fn new() -> Self {
        Self{start_time : Instant::now(),
            best_yet: 0.0, tresh : vec![], treshmax : vec![], timeout : 150.0, sliding : 0, registerName : String::new(), pruned : 0, spared : 0 }
    }

    pub fn lazy_nmcs_v3(&mut self, mut st: State, n : i8, r : f64, p : i8, p_budget : usize, heuristic_w : f64, verbose : bool) -> State{
        //println!("nmcs level {}", n);
        let mut best_state: State = st.clone(); //State::new();
        let mut best_state_score = -1.0; //best_state.score();
        while !st.terminal(){
            let moves = st.legal_moves();
            //println!("{}", self.start_time.elapsed().as_secs_f64());
            //println!("number of moves : {}", moves.len());
            if moves.len() == 0 || (self.start_time.elapsed().as_secs_f64() > self.timeout && self.timeout > 0.0) {return best_state;}

            let mut candidate_states : Vec<WS> = Vec::new();
            //println!("eval");

            let mut width= p_budget;
            if p_budget == 0 {
                width = moves.len();
            }

            for i in 0..width {
                //println!("{}", i);
            //for &mv in &moves{
                let mut index = ((moves.len() as f64)*rand::random::<f64>()) as usize;
                if p_budget == 0 {
                    index = i;
                }
                let mv = moves[index];


                if self.start_time.elapsed().as_secs_f64() > self.timeout && self.timeout > 0.0 {return best_state;}

                let mut new_st = st.clone();
                new_st.play(mv);

                if n <= 1 {
                    new_st = playout(new_st, heuristic_w);

                    candidate_states.push(WS{w : 0.0, s : new_st});
                }else{
                    let mut eval_Scores = 0.0;
                    let mut evalScore = 0.0;
                    if p != 0 {

                        /*
                        for _ in 0..p{
                            let mut eval_st = new_st.clone();
                            let handle = thread::spawn( move || {
                                eval_st = playout(eval_st, heuristic_w);
                                WS{w : eval_st.score(), s : eval_st}
                            });
                            handles.push(handle);
                        }

                        for h in handles {
                            let res = h.join().unwrap();
                            eval_Scores += res.w;


                            if res.w > best_state_score {
                                best_state = res.s;
                                best_state_score = res.w;
                                if best_state_score > self.best_yet {
                                    println!("in eval");
                                    println!(" pruned : {}  spared : {}", self.pruned, self.spared);
                                    self.best_yet = best_state_score;
                                    let elapsed = self.start_time.elapsed().as_secs_f64();
                                    writeLine(elapsed.to_string() + " " + &*best_state_score.to_string()+ "\n", self.registerName.clone());
                                    if verbose {println!("best score yet : {} after {}", best_state_score, elapsed);}


                                }
                            }
                        }
                         */
                        for _ in 0..p{
                            let mut eval_st = new_st.clone();
                            eval_st = playout(eval_st, heuristic_w);

                            //println!("leng {}", eval_st.seq.len());

                            let eval_sc = eval_st.score();
                            eval_Scores += eval_sc;


                            if eval_sc > best_state_score {
                                best_state = eval_st;
                                best_state_score = eval_sc;
                                if best_state_score > self.best_yet {
                                    if verbose {
                                        println!("in eval");
                                        println!(" pruned : {}  spared : {}", self.pruned, self.spared);
                                    }
                                    self.best_yet = best_state_score;
                                    let elapsed = self.start_time.elapsed().as_secs_f64();
                                    writeLine(elapsed.to_string() + " " + &*best_state_score.to_string()+ "\n", self.registerName.clone());
                                    if verbose {println!("lazyNMCSv3 best score yet : {} after {}", best_state_score, elapsed);}

                                    /*
                                    if best_state.reached_best_score {
                                        //on arrête tout
                                        if verbose{ println!("reached best score !!!");}
                                        return best_state
                                    }

                                     */



                                }
                            }
                        }

                        evalScore = eval_Scores/(p as f64);
                    }else{
                        evalScore = new_st.score();
                    }

                    if self.treshmax.len() < st.seq.len()+1 {
                        self.treshmax.push(evalScore);
                    }else{
                        if self.treshmax[st.seq.len()] < evalScore {
                            self.treshmax[st.seq.len()] = evalScore;
                        }
                    }

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

                    candidate_states.push(WS{w : evalScore, s : new_st});

                }
            }
            //println!("fin eval");
            //println!("{}", candidate_states.len());
            //println!("{}", n);
            for candidate in candidate_states {
                let mut new_st = candidate.s;
                if n > 1 {
                    //println!(" score {}", candidate.w);
                    //println!(" tresh {}", self.tresh[st.seq.len()].val);
                    if candidate.w < self.tresh[st.seq.len()].val + (self.treshmax[st.seq.len()]-self.tresh[st.seq.len()].val) * r {
                        //println!("prunning! {}", self.pruned);
                        self.pruned += 1;
                        new_st = playout(new_st, heuristic_w);
                    }else{
                        //println!("spared! {}", self.spared);
                        self.spared += 1;
                        new_st = self.lazy_nmcs_v3(new_st, n-1, r, p, p_budget, heuristic_w, verbose);
                    }
                }


                let new_st_score = new_st.score();


                if new_st_score > best_state_score {
                    best_state = new_st;
                    best_state_score = new_st_score;
                    if best_state_score > self.best_yet {
                        if verbose {println!(" pruned : {}  spared : {}", self.pruned, self.spared);}
                        self.best_yet = best_state_score;
                        let elapsed = self.start_time.elapsed().as_secs_f64();
                        writeLine(elapsed.to_string() + " " + &*best_state_score.to_string()+ "\n", self.registerName.clone());
                        if verbose {println!("lazyNMCSv3 best score yet : {} after {}", best_state_score, elapsed);}

                        /*
                        if best_state.reached_best_score {
                            //on arrête tout
                            if verbose{ println!("reached best score !!!");}
                            return best_state
                        }

                         */


                    }
                }
                /*
                if best_state.reached_best_score {
                    //on arrête tout
                    if verbose{ println!("reached best score !!!");}
                    return best_state
                }

                 */
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

pub fn launch_lazy_nmcs_v3(init_stat : State, level : i8, r :f64, p : i8, p_budget : usize,  heuristic_w : f64, timeout : f64, sliding : i32, verbose : bool, registerName : String) -> State{
    let mut expe = lazyNMCSv3::new();
    expe.timeout = timeout;
    expe.sliding = sliding;
    expe.registerName = registerName;
    //let mut st = State::new();

    let mut st = expe.lazy_nmcs_v3(init_stat, level, r, p, p_budget,  heuristic_w, verbose);

    return st;
}
