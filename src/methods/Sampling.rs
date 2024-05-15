use crate::models::SMILESgen::{State, Move};

use std::collections::HashMap;
use rand::random;
use crate::tools::calc::softmaxChoice;
use std::time::Instant;
use nalgebra::max;
use crate::tools::resultSaver::writeLine;

pub fn playout(mut st: State, heuristic_w : f64) -> State{

    let mut best_state: State = st.clone();
    let mut best_state_score = best_state.score();

    while !st.terminal() {
        //println!(" {:?}", st.coalitions);
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

pub fn Sample(init_state : State,  maxExp : i32, heuristic_w : f64, timeout : f64, registerName : String) -> State{
    //let mut inist = State::new();
    let mut best_Score = 0.0;
    let mut best_state = init_state.clone();
    let expeTODO = if maxExp == 0 {1} else {maxExp};
    let mut i = 0;
    let mut start_time = Instant::now();

    while i < expeTODO && !(start_time.elapsed().as_secs_f64() > timeout && timeout > 0.0) {
        let mut new_st = playout(init_state.clone(), heuristic_w);
        if new_st.score() > best_Score {
            best_state = new_st.clone();
            best_Score = new_st.score();
            writeLine(start_time.elapsed().as_secs_f64().to_string() + " " + &*best_Score.to_string()+ "\n", registerName.clone());

        }

        //HP model et SMILEGEN
        if new_st.reached_best_score {
            //on arrête tout
            //if verbose{ println!("reached best score !!!");}
            //println!("reached best score !!!");
            return new_st
        }
    }

    return best_state;
}