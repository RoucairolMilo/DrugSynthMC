use nalgebra::min;
use crate::methods::BFS;
//use crate::models::conjectures::GenerateGraph::{State, Move};
//use crate::models::Network::{State, Move};
//use crate::models::conjectures::conjecture2p4::{State, Move};
use crate::models::conjectures::extremal::{State};
//use crate::models::nonoGramSolver::{State};
use crate::tools::calc::softmaxChoice;
use crate::methods::NMCS;
use std::time::Instant;
use crate::tools::resultSaver::writeLine;
use std::collections::HashMap;

static SKIP_REPEATING_SCORES : bool = false; // garder à false, ne marche aps sur le 289 par exemple

#[derive(Clone)]
pub struct WS{ //Weight-State
    pub w : f64,
    pub s : State

}


pub fn insertDicho(l : &Vec<WS>, node : &WS) -> usize{
    let mut i = l.len()/2;
    let mut mi = 0;
    let mut ma = l.len()-1;
    //println!("entrée insert");
    while (i!=0 && l[i-1].w > node.w) || l[i].w < node.w{
        if l[i].w == node.w {
            return i;
        }
        if l[i].w < node.w{
            mi = i + 1;
            if mi> ma{
                return mi;
            }
        }else{
            ma = i;
        }
        i = (mi as f64/2.0 + (ma as f64)/2.0) as usize
    }
    //println!("sortie insert");
    return i;
}

pub fn playout(mut st: State, heuristic_w : f64) -> State{

    //super greedy pour le HP model
    /*
    while !st.terminal() {

        let moves = st.iterative_depth_greedy_moves();
        //println!("{}", moves.len());
        if moves.len() == 0 {return st;}
        for mv in moves {
            //println!("mv : {}", mv.dir );
            st.play(mv);
            //println!("sc : {}", st.score() );
        }
    }

     */


    //println!("entrée playout");
    /*
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
    }

     */

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

    if State::CONSIDER_NON_TERM{
        return best_state;
    }
    return st;

    //println!("sortie playout");
    return st;
}

pub fn BFS(inist : State,heuristic_w : f64, p:i32, timeout : f64, registerName : String) -> State{

    if p< 0 {
        //println!("attention, utilisation des scores des état non finaux au lieu des scores de playouts pour déterminer la valeur d'un noeud")
    }

    let mut st = inist.clone();
    let mut start_time = Instant::now();

    let mut open_nodes = Vec::new();
    open_nodes.push(WS{w : 0.0, s : st.clone()});

    let mut best_score_yet = 0.0;
    let mut best_state_yet = st.clone();
    let mut opened_nodes = 0;

    let mut visitedScores: HashMap<i64, bool> = HashMap::new();

    while open_nodes.len() != 0 {
        //println!("open nodes : {}", open_nodes.len());


        opened_nodes += 1;
        //println!("{}", start_time.elapsed().as_secs_f64());
        if start_time.elapsed().as_secs_f64() > timeout && timeout > 0.0 {return best_state_yet;}

        let mut node = open_nodes.pop().unwrap();
        //
        //println!("node opened with score : {}", node.w);
        //println!("{}", node.s.adj_mat);

        if SKIP_REPEATING_SCORES {
            while visitedScores.contains_key(&((node.w*10000000000.0) as i64)) && node.w != 0.0 && open_nodes.len() != 0 {
                node = open_nodes.pop().unwrap();
            }
            visitedScores.insert(((node.w*10000000000.0) as i64), true);
        }


        for m in node.s.legal_moves() {
            let mut new_state = node.s.clone();
            new_state.play(m);

            if p >= 0 {

                //let mut expe = NMCS::NMCS::new();
                //let mut best_playout_state = expe.nmcs(new_state.clone(), 1, heuristic_w, false);

                let mut best_playout_state = playout(new_state.clone(), heuristic_w);
                //println!("playout !");
                //while true {}

                let mut best_playout_state_score = best_playout_state.score();
                for i in 0..p{

                    //let mut expe = NMCS::NMCS::new();
                    //let mut playout_state = expe.nmcs(new_state.clone(), 1, heuristic_w, false);

                    let mut playout_state = playout(new_state.clone(), heuristic_w);

                    let playout_state_score = playout_state.score();

                    if playout_state_score > best_playout_state_score {
                        best_playout_state = playout_state.clone();
                        best_playout_state_score = playout_state_score;
                    }
                }
                if best_playout_state_score > best_score_yet{
                    best_score_yet = best_playout_state_score;
                    best_state_yet = best_playout_state.clone();
                    writeLine(start_time.elapsed().as_secs_f64().to_string() + " " + &*best_score_yet.to_string()+ "\n", registerName.clone());
                    //println!("BFS record battu ! {}", best_score_yet);
                    println!("BFS record battu : {} after {}", best_score_yet, start_time.elapsed().as_secs_f64());

                }

                let new_ws = WS{w : best_playout_state_score, s : new_state};
                let mut i = 0;
                if open_nodes.len() != 0 {
                    i = insertDicho(&open_nodes, &new_ws);
                }
                open_nodes.insert(i, new_ws);
            }else{
                let smoothed_score = new_state.smoothedScore(); //une estimation de lintêret de continuer sur cet enfant
                if smoothed_score > best_score_yet{
                    let new_state_score = new_state.score();
                    if new_state_score > best_score_yet{
                        best_score_yet = new_state_score;
                        best_state_yet = new_state.clone();
                        writeLine(start_time.elapsed().as_secs_f64().to_string() + " " + &*best_score_yet.to_string()+ "\n", registerName.clone());
                        //println!("BFS record battu ! {}", best_score_yet);
                        println!("BFS record battu : {} after {}", best_score_yet, start_time.elapsed().as_secs_f64());
                        //println!("{}", best_state_yet.adj_mat);
                        //for i in 0..open_nodes.len() {print!("{} ", open_nodes[open_nodes.len() -1 - i].w);}
                        /*
                        if best_state_yet.reached_best_score {
                            return best_state_yet;
                        }

                         */

                        //println!(" ");
                    }
                }
                if !new_state.terminal() {
                    //println!("{}, {}", m.from, m.to);
                    let new_ws = WS{w : smoothed_score, s : new_state};

                    let mut i = 0;
                    if open_nodes.len() != 0 {
                        i = insertDicho(&open_nodes, &new_ws);
                    }
                    open_nodes.insert(i, new_ws);
                }

            }
        }
        println!("open node {}", open_nodes.len());
        println!("open node {}", best_state_yet.terminal());
    }
    //println!(" exhausted all nodes : {}", opened_nodes);
    return best_state_yet;
}


pub fn launch_bfs(init_stat : State, heuristic_w : f64, p:i32, timeout : f64, registerName : String) -> State{
    return BFS(init_stat, heuristic_w, p, timeout, registerName);
}
