use crate::models::SMILESgen::{State, Move};

use std::collections::HashMap;
use rand::random;
use crate::tools::calc::softmaxChoice;
use std::time::Instant;
use nalgebra::max;
use crate::tools::resultSaver::writeLine;

#[derive(Clone)]
struct aver{ //WeightMove
pub n : i32,
    pub a : f64
}

pub fn playout(mut st: State, heuristic_w : f64) -> State{

    let mut best_state: State = st.clone();

    let mut best_state_score = 0.0;
    if State::CONSIDER_NON_TERM || st.terminal() {
        best_state_score = best_state.score();
    }

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

pub fn UCT(inist : State, Cp : f64, maxExp :i32, heuristic_w :f64, verbose : bool, timeout : f64, registerName : String) -> State{
    let mut visits : HashMap<Vec<Move>, aver> = HashMap::new();
    let mut explo = 0;
    let mut bestScore = 0.0;
    let mut bestState = inist.clone();
    let mut start_time = Instant::now();


    while explo != maxExp {
        explo +=1;
        //println!("visits : {}", visits.len());
        //recharche à partir du premier noeud
        let mut st = inist.clone();
        let mut moveSeq: Vec<Move> = vec![];
        let mut leaf = false;

        while !leaf {
            if start_time.elapsed().as_secs_f64() > timeout && timeout > 0.0 {return bestState;}
            leaf = true;
            let mut best_UCT_score = 0.0;

            if !visits.contains_key(&moveSeq) {
                //println!("explore new");
                //jamais visité
                let mut new_st = playout(st.clone(), heuristic_w);
                let sc = new_st.score();
                if sc > bestScore{
                    bestScore = sc;
                    bestState = new_st.clone();

                    writeLine(start_time.elapsed().as_secs_f64().to_string() + " " + &*bestScore.to_string()+ "\n", registerName.clone());
                    //println!("UCT new best score A : {} in {}", bestScore, start_time.elapsed().as_secs_f64());

                    //HP model et SMILEGEN
                    if new_st.reached_best_score {
                        //on arrête tout
                        //if verbose{ println!("reached best score !!!");}
                        return new_st
                    }
                }

                if moveSeq.len() == 0 {
                    visits.insert(moveSeq.clone(), aver{n : 1, a : bestState.score()});
                }else{
                    for i in 1..moveSeq.len()+1 {
                        if visits.contains_key(&moveSeq[0..i]) {
                            let na = visits.get(&moveSeq[0..i]).unwrap();
                            visits.insert(Vec::from(&moveSeq[0..i]), aver{n : na.n+1, a : (na.a * na.n as f64 + sc)/(na.n as f64 + 1.0 )});
                        }else{
                            visits.insert(Vec::from(&moveSeq[0..i]), aver{n : 1, a : sc});
                        }
                    }
                }
            }else
            {

                let mut next_moveSeq = moveSeq.clone();
                let nombre_visite_pere = visits.get(&moveSeq).unwrap().n;
                let moves = st.legal_moves();

                if moves.len() == 0 || st.terminal() {
                    //feuille, remonte le resultat
                    //println!("exploit old");
                    //println!("{:?}, {}",st.SMILE, st.clone().score());
                    leaf = true;
                    let sc = st.score();
                    if sc > bestScore {
                        bestScore = sc;
                        bestState = st.clone();
                        writeLine(start_time.elapsed().as_secs_f64().to_string() + " " + &*bestScore.to_string()+ "\n", registerName.clone());
                        //println!("UCT new best score B : {} in {}", bestScore, start_time.elapsed().as_secs_f64());

                        //HP model et SMILEGEN
                        if st.reached_best_score {
                            //on arrête tout
                            //if verbose{ println!("reached best score !!!");}
                            return st
                        }
                    }
                    if moveSeq.len() == 0 {
                        visits.insert(moveSeq.clone(), aver{n : 1, a : bestState.score()});
                    }else{
                        for i in 1..moveSeq.len()+1 {
                            if visits.contains_key(&moveSeq[0..i]) {
                                let na = visits.get(&moveSeq[0..i]).unwrap();
                                visits.insert(Vec::from(&moveSeq[0..i]), aver{n : na.n+1, a : (na.a * na.n as f64 + sc)/(na.n as f64 + 1.0 )});
                            }else{
                                visits.insert(Vec::from(&moveSeq[0..i]), aver{n : 1, a : sc});
                            }
                        }
                    }


                }else{
                    //println!("progress down");
                    let mut best_move: Move = moves[0];

                    for mv in moves {
                        let mut nombre_visite_fils =0;
                        let mut ave_fils = 0.0;
                        let mut new_moveSeq = moveSeq.clone();
                        new_moveSeq.push(mv);

                        if visits.contains_key(&new_moveSeq) {
                            let temp = visits.get(&new_moveSeq).unwrap();
                            nombre_visite_fils = temp.n;
                            ave_fils = temp.a;
                        }

                        //choisis le noeud à la valeur la plus élevée
                        let mut UCT_score = f64::NEG_INFINITY;

                        if nombre_visite_fils != 0 {
                            UCT_score = ave_fils/bestScore + Cp * (2.0* (nombre_visite_pere as f64).ln() / (nombre_visite_fils as f64) ).sqrt();
                            //println!("explique {} et {} mv {} bs {}", ave_fils/bestScore, Cp * (2.0* (nombre_visite_pere as f64).ln() / (nombre_visite_fils as f64) ).sqrt(), mv.atom, bestScore);
                            //println!("{:?}",st.SMILE);
                        }else{
                            UCT_score = f64::INFINITY;
                        }
                        if best_UCT_score < UCT_score || (best_UCT_score == UCT_score && rand::random::<f64>() < 0.5 ) {

                            next_moveSeq = new_moveSeq.clone();
                            best_UCT_score = UCT_score;
                            best_move = mv.clone();
                            //println!("chose : {}", mv.atom);
                        }
                    }

                    // joue le meilleur coup
                    moveSeq = next_moveSeq;
                    st.play(best_move);
                    leaf = false;
                    if best_UCT_score != 0.0 {
                    }
                }
            }
        }
    }
    return bestState;
}


pub fn launch_UCT(init_state : State, Cp : f64, maxExp : i32, heuristic_w : f64,timeout : f64, registerName : String) -> State{
    //let mut inist = State::new();
    return UCT(init_state, Cp, maxExp, heuristic_w, true, timeout, registerName);
}
