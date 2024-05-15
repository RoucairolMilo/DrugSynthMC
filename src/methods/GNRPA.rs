//use crate::models::coalitionsB::{State, Move};
//use crate::models::Network::{State, Move};
use crate::models::Network::{State, Move};
use std::collections::HashMap;
use std::ptr::hash;
use std::time::Instant;
use crate::tools::resultSaver::writeLine;

//utiliser plutôt que NRPA
const T : f64 = 1.0; //température du GNRPA
const alpha : f64 = 1.0; //taux d'apprentissage

pub fn random_move(st : &State, moves : Vec<Move>, policy: &mut HashMap<Move, f64>, heuristic_w: f64) -> Move{
    let mut sum : f64 = 0.0;

    let mut hash_heuri: HashMap<Move, f64> = HashMap::new();

    for &mv  in &moves {

        if heuristic_w > 0.0{
            hash_heuri.insert(mv, st.heuristic(mv)* heuristic_w);
        }else{
            hash_heuri.insert(mv, 0.0);
        }

        match policy.get(&mv){
            Some(v) => sum += (v/T + hash_heuri.get(&mv).unwrap()).exp(),
            None => {
                //policy.insert(mv, *hash_heuri.get(&mv).unwrap());
                policy.insert(mv, 0.0);
                sum+= hash_heuri.get(&mv).unwrap().exp()
            }
        };
    }
    let stop = sum*rand::random::<f64>();
    sum = 0.0;
    for &mv in &moves {
        sum += (policy.get(&mv).unwrap()/T + hash_heuri.get(&mv).unwrap()).exp() ;
        if sum > stop {
            return mv;
        }
    }
    return moves[0];
}

pub fn adapt(mut policy: HashMap<Move, f64>, st : &mut State, ini_st : &State, heuristic_w : f64) -> HashMap<Move, f64>{
    let mut s = ini_st.clone();
    let mut polp: HashMap<Move, f64> = policy.clone();

    for best in & st.seq[..] {
        let moves = s.legal_moves();
        let mut z = 0.0;
        let mut hash_heuri: HashMap<Move, f64> = HashMap::new();
        for &m in &moves {
            if heuristic_w > 0.0{
                hash_heuri.insert(m, s.heuristic(m)* heuristic_w);
            }else{
                hash_heuri.insert(m, 0.0);
            }

            match policy.get(&m){
                Some(v) => z += (v/T + hash_heuri.get(&m).unwrap()).exp(),
                None => {
                    //policy.insert(m, hash_heuri.get(&m).unwrap().exp());
                    //policy.insert(m, 0.0);
                    //z += policy.get(&m).unwrap()}
                    z+= hash_heuri.get(&m).unwrap().exp()}
            };
        }

        for &m in &moves {
            let mut delta = 0.0;
            if &m == best {delta = 1.0}
            match polp.get(&m){
                Some(v) => polp.insert(m, v - alpha/T *  ((v/T + hash_heuri.get(&m).unwrap()).exp()/z - delta)),
                None => polp.insert(m, -alpha/T * (hash_heuri.get(&m).unwrap().exp()/z-delta))
            };
        }

        s.play(*best);
    }
    return polp;
}


pub struct GNRPA{
    pub best_yet : f64,
    pub timeout : f64,
    pub registerName : String,
    pub start_time : Instant,
}

impl GNRPA{
    pub fn new() -> Self{
        Self{ start_time : Instant::now(),best_yet: 0.0, timeout : -1.0, registerName : String::new() }
    }

    pub fn playout(&mut self, mut st : State, mut policy : HashMap<Move, f64>, heuristic_w: f64) -> State{
        let mut best_state: State = st.clone(); //State::new();
        let mut best_state_score = best_state.score();

        while !st.terminal() {
            let  moves: Vec<Move> = st.legal_moves();
            if moves.len() == 0 {
                return st;
            }

            let mv = random_move(&st, moves, &mut policy, heuristic_w);

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

    pub fn gnrpa(&mut self, mut ini_st: State, level : i8, mut policy: HashMap<Move, f64>, heuristic_w: f64, initial : bool) -> State {
        let mut st: State = ini_st.clone();
        let mut stscore : f64 = st.score();

        if level == 0 {

            let osef =  self.playout(st, policy, heuristic_w);

            return osef;
        }

        for i in 0..100 {
            if self.start_time.elapsed().as_secs_f64() > self.timeout && self.timeout > 0.0 {return st;}
            if initial {println!("GNRPA loop {}, best score : {} ", i, stscore);}
            let pol : HashMap<Move, f64> = policy.clone();
            let mut s = self.gnrpa( ini_st.clone(),level-1, pol, heuristic_w, false);

            let s_score = s.score();

            if stscore < s_score {
                st = s;
                stscore = s_score;

                if(self.best_yet < stscore){
                    let elapsed = self.start_time.elapsed().as_secs_f64();
                    self.best_yet = stscore;
                    println!("best score yet : {} after {}", stscore, elapsed);
                    writeLine(elapsed.to_string() + " " + &*stscore.to_string()+ "\n", self.registerName.clone());
                }
            }

            policy = adapt(policy, &mut st, &ini_st, heuristic_w);

            if level == 1 {
                //println!("best yet : {}, {}, {}",stscore, st.n_sommet, st.n_arete );
            }
        }

        return st;
    }
}










pub fn launch_gnrpa(ini_st : State, level : i8, heuristic_w: f64, timeout : f64, registerName : String) -> State {
    let mut expe = GNRPA::new();
    expe.timeout = timeout;
    expe.registerName = registerName;

    let mut policy = HashMap::new();
    let st = expe.gnrpa(ini_st.clone(), level, policy, heuristic_w, true);
    return st;
}