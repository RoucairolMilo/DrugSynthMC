use crate::models::HPmodel::{State, Move};
use crate::tools::calc::softmaxChoice;

impl State {
    pub fn noHeuristic(st : State) -> f64{
        return 1.0;
    }
}

pub struct PermParam{
    pub maxSamples : i32, //combien de playouts
    pub k : i32, //enrichir/cloner en cb d'enfants
    pub lowW : f64, //seuil de prunning
    pub highW : f64, //seuil d'enrichissement
    pub boltzFact : f64, //1.3 est bien d'après l'article
}


impl PermParam{
    pub fn new() -> Self{
        Self{maxSamples : 100000, k : 1, lowW : 0.2, highW : 1.0, boltzFact : 1.0}
    }
}

#[derive(Clone)]
struct WSP{ //Weight-State-Proba
    pub w : f64,
    pub s : State,
    pub p : f64
}

#[derive(Clone)]
struct WM{ //WeightMove
    pub w : f64,
    pub m : Move
}

pub fn PERM(mut st : State, p : PermParam) -> State{

    let mut best = 0.0;
    let mut bestState= st.clone();
    let mut samples = 0;
    let mut candidates: Vec<WM>= Vec::new();
    let mut candidatesW : Vec<f64>= Vec::new();
    let mut stack : Vec<WSP> = Vec::new();

    let mut Z = 1.0;
    let mut all_trials_w = Vec::new();
    for _ in 0..st.molecule.len()+1{
        all_trials_w.push(vec![]);
    }

    stack.push(WSP{w : 1.0, s : st.clone(), p : 1.0});
    stack.push(WSP{w : 1.0, s : st, p : 1.0});

    while stack.len() > 0 && samples < p.maxSamples {
        samples += 1;
        let mut ws: WSP = stack[0].clone();
        stack.remove(0);

        //for t in &stack{print!("{} ",t.s.score());}
        //println!("");
        //println!("{} ",stack.len());

        while !ws.s.terminal() {
            candidates = Vec::new();
            candidatesW = Vec::new();
            let mut sumexp = 0.0;
            for m in ws.s.legal_moves() {
                let w = ws.s.heuristic(m)*p.boltzFact;
                candidates.push(WM{ w : w, m : m});
                candidatesW.push(w);
                sumexp += w.exp();
            }

            if candidates.len() == 0 {
                //bloqué, on ouvre une autre copie
                break;
            }
            else{
                let ind = softmaxChoice(candidatesW.clone());
                let mv  = candidates[ind].m;
                ws.s.play(mv);
                ws.w *= sumexp;
                ws.p *= candidates[ind].w.exp()/sumexp;

                if all_trials_w[ws.s.seq.len()].len() != 0 {
                    Z = all_trials_w[ws.s.seq.len()].iter().sum::<f64>() / (all_trials_w[ws.s.seq.len()].len() as f64);
                }else{
                    Z = 0.0;
                }
                all_trials_w[ws.s.seq.len()].push((p.boltzFact*ws.s.score()).exp()/ws.p);

                if ws.w > p.highW*Z && Z!=0.0{
                    ws.w *= 1.0/(1.0+(p.k as f64));
                    for i in 0..p.k{
                        stack.insert(0, WSP{w : ws.w, s : ws.s.clone(), p : ws.p});
                    }
                }

                if ws.w < p.lowW*Z && Z != 0.0{
                    if rand::random::<f64>() < 0.5 {
                        //pruning
                        break;
                    }else{
                        ws.w *= 2.0;
                    }
                }
            }
        }

        if candidates.len() == 0{
            //println!("pas fini");
            //println!("stack length : {}", stack.len());
        }else{
            let sc = ws.s.score();
            if sc > best{
                best = sc;
                bestState = ws.s.clone();
                println!("best yet : {}", sc);
                //println!("stack length : {}", stack.len());
            }
        }
    }
    return bestState;
}

/*
fn main() {
    println!("Hello, world!");

    let mut p = PermParam::new();
    p.lowW = 0.5;
    p.highW = 5.0;
    let mut st = HPmodel::State::new();
    //st.molecule = "HHHHHHHH".parse().unwrap();
    PERM::PERM(st, p);
}
 */