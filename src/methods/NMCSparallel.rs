use crate::models::SMILESgen::{State, Move};
use std::thread;
use std::thread::JoinHandle;
use crate::tools::calc::softmaxChoice;

pub struct standsc{
    pub sc : f64,
    pub st : State
}



pub fn playout(mut st: State, heuristic_w : f64) -> State{
    while !st.terminal() {
            let moves = st.legal_moves();
            if moves.len() == 0 { return st; }
            let mut i = ((moves.len() as f64)*rand::random::<f64>()) as usize;
            if heuristic_w != 0.0 {
                let mut weights = Vec::new();
                for &m in &moves {
                    weights.push(heuristic_w * st.heuristic(m));
                }
                i = softmaxChoice(weights);
            }
            let mv = moves[i as usize];
            st.play(mv);
        }
        return st;
    }

    pub fn nmcs(mut st: State, n : i8, heuristic_w : f64) -> State{
        let mut best_state: State = State::new();
        let mut best_state_score = best_state.score();

        while !st.terminal(){
            let moves = st.legal_moves();
            if moves.len() == 0 { return st; }
            let mut handles : Vec<JoinHandle<standsc>> = Vec::new();

            let i = 0;
            for &mv in &moves{
                let mut new_st = st.clone();

                let handle = thread::spawn( move || {
                    new_st.play(mv);
                    if n == 1 {
                        new_st = playout(new_st, heuristic_w);
                    }else{
                        new_st = nmcs(new_st, n-1, heuristic_w);
                    }
                    //println!("s : {}",newSt.score() );
                    let new_st_score = new_st.score();

                    standsc{ sc: new_st_score, st : new_st }
                });
                handles.push(handle);
            }

            for h in handles {
                let res = h.join().unwrap();
                if res.sc > best_state_score {
                    best_state = res.st;
                    best_state_score = res.sc;
                }
            }
            //println!("a : {}",st.seq.len() );
            //println!("b : {}",bestState.seq.len() );
            //println!("best score : {}, {}, {}", best_state_score, best_state.n_sommet, best_state.n_arete );
            st.play(best_state.seq[st.seq.len()]);
        }
        return st;
    }


pub fn launch_para_nmcs(level : i8, heuristic_w : f64) -> State{

    let mut st = State::new();
    let mut bestScoreYet = standsc{sc : 0.0, st : State::new() };
    st = nmcs(st, level, heuristic_w);
    return st;
}
