

extern crate nalgebra;
use std::collections::{HashMap, HashSet};
use nalgebra::{ArrayStorage, Const, Matrix, max, Vector3};

use serde::Deserialize;
use serde_json::Result;
use std::fs;
use std::clone::Clone;

use crate::tools::NNreader::{Prediction, Model};
//use crate::tools::pythonBind::{NNholder};


//S(=O)(=O) -> U
//S(=O) -> M
//Cl -> L
//C(F)(F)(F) -> W
const ATOMS: [char; 10] = ['C', 'O', 'N', 'F', 'S', 'P', 'U', 'M', 'L', 'W'];
const NUMBERS: [char; 9] = ['1', '2', '3', '4', '5', '6', '7', '8', '9'];
const NGRAM_LEN: usize = 3; //at least 1
const PRUNING_MOVE_TR: f64 = 0.001; //frequency under which we prune the move, deemed as too rare, 0.01 for ngrams and neural, 0.001 with the sulfur now? at least 0.0 when using ngrams to remove the moves with a score of 0.0
const MAX_RING_LEN: usize = 7;
const MIN_RING_LEN: usize = 5;
const HEURISTIC_MODE: &str = "ngram"; // "ngram" "neural"


#[derive(Debug, Deserialize)]
struct NA {
    #[serde(flatten)]
    nextAtom: HashMap<char, f64>,
}

lazy_static! {
    static ref legalBonds: HashSet<(char, char, i32)> = {
        let mut m = HashSet::new();
        m.insert(('C', 'C', 1));
        m.insert(('C', 'O', 1));
        m.insert(('C', 'N', 1));
        m.insert(('C', 'F', 1));
        m.insert(('C', 'S', 1));
        m.insert(('C', 'P', 1));
        m.insert(('C', 'B', 1));
        m.insert(('C', 'c', 1));
        m.insert(('C', 'I', 1));
        m.insert(('O', 'N', 1));
        m.insert(('N', 'N', 1));
        m.insert(('N', 'S', 1));
        m.insert(('N', 'P', 1));
        m.insert(('P', 'S', 1));

        m.insert(('C', 'C', 2));
        m.insert(('C', 'O', 2));
        m.insert(('C', 'N', 2));
        m.insert(('C', 'S', 2));
        m.insert(('O', 'N', 2));
        m.insert(('O', 'S', 2));
        m.insert(('N', 'N', 2));
        m.insert(('N', 'P', 2));
        m.insert(('S', 'P', 2));

        m
    };


    static ref prior : Model = {
        let m = Model::from_dir("Neural/SMILESexplicit_shortcuts");
        //let m = Model::from_dir("Neural/testActix");

        let mo = m.unwrap();
        //println!("yabadabadooo {:?}", mo.predict(vec!["C", "O", "C"]).unwrap());

        mo
    };



    /*
    static ref prior : NNholder = {
        let m = pythonBind::NNholder::from_dir("Neural/SMILES");


        println!("yabadabadooo {:?}", m.predict(vec!["C", "O", "C"]).unwrap());

        m
    };

     */

    static ref ngrams : HashMap<String, NA> = {
        //let path = "ngrams/ngrams.json";
        let path = "ngrams/fda_ngrams_shortcuts_cycles.json";
        //let path = "ngrams/zinc_ngrams_1.json";
        let ngram_string = fs::read_to_string(path).expect("Unable to read file");
        //println!("{}",ngram_string);
        let ngrams_data: HashMap<String, NA> = serde_json::from_str(&ngram_string).unwrap();

        ngrams_data
    };








}



#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub struct Move{
    pub atom : char,
    pub doubleLink : bool,
    pub nesting : bool,
    pub closeNesting : bool,
    pub cycle : usize
}

#[derive(Clone)]
pub struct State{
    pub SMILE : Vec<char>,
    pub nestingOpenAtom : Vec<char>,
    pub nestingOpenCovalence : Vec<i32>,
    pub nestingCycleToClose : Vec<usize>,
    pub openCycles : Vec<usize>,
    pub openCyCosts : Vec<usize>,
    pub closedCycles : usize,
    pub seq : Vec<Move>,
    pub open_nesting_ASAP : bool,
    pub finish_ASAP : bool,
    pub reached_best_score :bool,

    pub target_OtoC_ratio : f64,
    pub target_NtoC_ratio : f64,

    pub storedPrior : Prediction
}

impl State{
    pub const CONSIDER_NON_TERM: bool = false;
    pub const BEST_POSSIBLE_SCORE: f64 = 1002.0;

    pub fn new() -> Self {

        Self{
            SMILE : vec!['C'],
            nestingOpenAtom : vec!['C'],
            nestingOpenCovalence : vec![4],
            nestingCycleToClose : vec![0],
            openCycles : Vec::new(),
            openCyCosts : Vec::new(),
            closedCycles : 0,

            seq : Vec::new(),
            open_nesting_ASAP : false,
            finish_ASAP : false,
            reached_best_score : false,

            target_OtoC_ratio : -1.0,
            target_NtoC_ratio : -1.0,

            storedPrior : Prediction{label : vec![], confidence: vec![] }
        }
    }



    pub fn play(&mut self, m : Move){

        self.storedPrior = Prediction{label : vec![], confidence: vec![] };

        let mut addition : char = ' ';
        let s = self.nestingOpenCovalence.len().clone();

        if m.doubleLink && m.cycle == 0 {
            addition = '=';


            //pour les double liaisons après une parenthèse : affecte le nesting d'avant, sinona ffecte le nesting courant
            if self.SMILE[self.SMILE.len()-1] == '('{
                self.nestingOpenCovalence[s -2] -= 1;
            }else{
                self.nestingOpenCovalence[s -1] -= 1;
            }
        }

        if m.nesting {
            addition = '(';
            self.nestingOpenCovalence[s -1] -= 1;
            self.nestingOpenCovalence.push(1);
            self.nestingOpenAtom.push('X');

            //self.nestingParentIndex.push(n);
            if self.open_nesting_ASAP {
                self.open_nesting_ASAP = false;
                //println!("{:?}", self.openCycles);
                //println!("{:?}", self.SMILE);
                self.nestingCycleToClose.push(self.openCycles[self.openCycles.len()-1])
            }else{
                self.nestingCycleToClose.push(0);
            }

        }

        if m.closeNesting {
            addition = ')';
            self.nestingOpenCovalence.pop();
            self.nestingOpenAtom.pop();
            self.nestingCycleToClose.pop();
        }

        if m.cycle != 0 {
            addition = char::from_digit(m.cycle as u32, 10).unwrap();

            if self.openCycles.contains(&m.cycle) {
                self.closedCycles += 1;
                self.openCycles.retain(|&x| x != m.cycle);

                for i in 0..self.nestingCycleToClose.len() {
                    if self.nestingCycleToClose[i] == m.cycle {
                        self.nestingCycleToClose[i] = 0;

                    }
                }

            }else{
                self.open_nesting_ASAP = true;
                self.openCycles.push(m.cycle);

            }

            self.nestingOpenCovalence[s -1] -= 1;
            if m.doubleLink {
                self.nestingOpenCovalence[s -1] -= 1;
            }
        }

        if m.atom != ' ' {
            addition = m.atom;

            let mut new_covalence = 0;
            if m.atom == 'C' { new_covalence = 3; } //toujours une liaison de moins
            if m.atom == 'N' { new_covalence = 2; }
            if m.atom == 'O' { new_covalence = 1; }
            if m.atom == 'F' { new_covalence = 0; }
            if m.atom == 'S' { new_covalence = 1; } //considéré comme un oxygene dans le cas S
            if m.atom == 'U' { new_covalence = 1; }
            if m.atom == 'M' { new_covalence = 1; }
            if m.atom == 'L' { new_covalence = 0; }
            if m.atom == 'W' { new_covalence = 0; }
            //if m.atom == 'P' { new_covalence = 4; }
            //println!("quoi ??? {} {}", m.atom, new_covalence);



            for i in 0..self.SMILE.len() {

                let last_char = self.SMILE[self.SMILE.len()-1-i];

                if last_char == '='{
                    //self.nestingOpenCovalence[s -1] -= 1;
                    new_covalence-=1;
                }else{
                    if last_char != '(' {
                        if last_char != ')'{
                            break
                        }
                    }else{
                        //self.nestingOpenCovalence.push(0);
                        //self.nestingParentIndex.push(0);
                    }

                }
            }

            self.nestingOpenCovalence[s -1] = new_covalence;
            self.nestingOpenAtom[s -1] = m.atom;

            //self.nestingParentIndex[s -1] = self.SMILE.len();

        }

        self.SMILE.push(addition);

        //if self.SMILE.len() > 10 { self.finish_ASAP = true }
        if self.SMILE.len() > 20 { self.finish_ASAP = true }

        //self.SMILE = *new_SMILE.clone();
        self.seq.push(m);
    }

    pub fn legal_moves(&mut self) -> Vec<Move>{

        let mut vec :Vec<Move> = Vec::new();

        let last_addition = self.SMILE[self.SMILE.len()-1];
        let mut last_last_addition = 'x';
        if self.SMILE.len() > 2{
            last_last_addition = self.SMILE[self.SMILE.len()-2];
        }
        //je regarde si il y a un niveau de nesting ouvert qui pourrait continuer après, s'il n'y en a pas, finir par un fluor, ou une double liason et un souffre ou oxygène alors qu'il y a un cycle d'ouvert est interdit
        let mut could_prevent_cycle_completion = true;

        for i in 0..self.nestingOpenCovalence.len()-1 {
            if self.nestingCycleToClose[self.nestingOpenCovalence.len()-1 - i] != 0 {
                break
            }
            if self.nestingOpenCovalence[self.nestingOpenCovalence.len()-2 - i] != 0 {could_prevent_cycle_completion = false }
        }

        //println!("could {}", could_prevent_cycle_completion);
        //println!("a {:?}", self.nestingOpenCovalence);
        //println!("b {:?}", self.nestingCycleToClose);
        let mut can_play_end_cycle = false;
        let mut must_close_cycle = false;

        if self.openCycles.len() >= 1 {
            let (cycle_len, ar, ri, proper_atoms , bc ) = Self::backtrackCycle(self.SMILE.clone(), char::from_digit(self.openCycles[self.openCycles.len()-1] as u32, 10).unwrap());
            let bond_cost = bc as i32;
            //on termine avec un cycle quand il y ne reste plus qu'une liaison seulement quand il ne reste plus dautre cycle à boucler
            if !self.finish_ASAP || (self.openCycles.len() == 1 || self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 2 + bond_cost-1 ) {
                //fermeture d'un cycle
                //et si le cycle se termine par un "=" ???
                if last_addition != '=' && !(last_addition == '(' && last_last_addition == '=' )  {
                    //interdit de fermer un cycle si ça bloque tout ensuite
                    if !(self.openCycles.len() > 1 && could_prevent_cycle_completion && self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] <= 1 + bond_cost-1) && self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 1 + bond_cost -1 {
                        /*
                        for i in &self.openCycles {
                            let cycle_len = Self::backtrackCycle(self.SMILE.clone(), char::from_digit(*i as u32, 10).unwrap()).0;
                            if cycle_len > 4 {
                                let mv = Move { atom: ' ', doubleLink : false, nesting: false, closeNesting: false, cycle: *i };
                                vec.push(mv);
                                can_play_end_cycle = true;
                            }
                            if cycle_len > 6 {
                                must_close_cycle = true;
                            }

                        }

                         */
                        if self.openCycles.len() > 0 {

                            if cycle_len >= MIN_RING_LEN && proper_atoms >= 4 {

                                let mut mv = Move { atom: ' ', doubleLink : false, nesting: false, closeNesting: false, cycle : self.openCycles[self.openCycles.len()-1] };
                                if bond_cost == 2 {
                                    mv.doubleLink = true;

                                }
                                vec.push(mv);
                                can_play_end_cycle = true;
                            }
                            if cycle_len >= MAX_RING_LEN {
                                must_close_cycle = true;
                            }
                        }
                    }
                }
            }

            //pour eviter els cycles un peu trop longs
            if must_close_cycle && can_play_end_cycle && (!could_prevent_cycle_completion || self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 2+bond_cost-1) && self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 1 + bond_cost-1  {
                return vec;
            }
        }


        //fermeture d'un nesting, toujours autorisé sauf après un nesting ou un double lien ou avant d'avoir refermé une boucle
        if self.nestingOpenCovalence.len() != 1 && self.nestingCycleToClose[self.nestingCycleToClose.len() -1 ] == 0 && !self.open_nesting_ASAP {
            if last_addition != '=' && last_addition != '(' && !(could_prevent_cycle_completion && self.openCycles.len() > 0) {
                let mv = Move { atom: ' ', doubleLink : false, nesting: false, closeNesting: true, cycle: 0 };
                vec.push(mv);
            }
        }

        if must_close_cycle && self.nestingCycleToClose[self.nestingCycleToClose.len() -1] != 0 {

            return vec;
        }

        //AIzynthfinder interdit les "=(", seulement les "(=" sont acceptés
        if self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 1 || last_addition == '(' {
            //ouverture d'un nesting, interdit en finish ASAP
            if  last_addition != '(' && last_addition != '=' &&  (!self.finish_ASAP || (must_close_cycle && !can_play_end_cycle)) {
                let mv = Move{atom : ' ', doubleLink : false, nesting : true, closeNesting : false, cycle : 0 };
                vec.push(mv);
            }
            if self.open_nesting_ASAP && ATOMS.contains(&last_addition)  {
                //rien, c'est foireux ça prend pas les nesting en compte
                let mv = Move{atom : ' ', doubleLink : false, nesting : true, closeNesting : false, cycle : 0 };
                return vec![mv];
            }
        }

        if self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 1 || last_addition == '(' {
            if !self.finish_ASAP || last_addition == '=' || last_addition == '(' || (self.openCycles.len() > 0) || (self.open_nesting_ASAP && !ATOMS.contains(&last_addition)) { //(self.openCycles.len() > 0 && !can_play_end_cycle) retiré la clause sur la fin de cycle pour éviter d'avoir trop de cycle de taille 5 forcés
                for i in ['C', 'O', 'F', 'N', 'S', 'U', 'M', 'L', 'W'] { //['C', 'O', 'N', 'F', 'P', 'S'] //['C', 'O', 'F', 'N', 'S'] trop de S
                    //println!("de ? {}", could_prevent_cycle_completion);

                    if i == 'O' {
                        //println!("condi A {}", self.finish_ASAP && last_addition != '=' && !could_prevent_cycle_completion);
                        //println!("condi B {}", self.openCycles.len() > 0 && could_prevent_cycle_completion);
                        //println!("condi C {}", self.finish_ASAP && last_addition == '=' && self.nestingCycleToClose[self.nestingCycleToClose.len() -1 ] != 0);

                        //println!("last {} nesting {:?} ASAP {}", last_addition, self.nestingCycleToClose, self.finish_ASAP);

                        //println!("condi D {}", last_addition == '=' && self.open_nesting_ASAP);
                    }

                    if !((i == 'F' || i == 'L' || i == 'W') && (self.nestingCycleToClose[self.nestingCycleToClose.len() -1 ] != 0 || (self.openCycles.len() > 0 && could_prevent_cycle_completion) || self.open_nesting_ASAP || last_addition == '=' || (last_addition == '(' && last_last_addition == '=' )) ) && //on saute aussi si c'est F et qu'on est au top niveau et qu'il reste des boucles à fermer
                        //on interdit S et O si on est en finish asap sauf si ça termine un = et que ça n'empêche aps de terminer les cycles, ou si pourrait bloquer la fin des cycles
                        !( (i == 'O' || i == 'U' || i == 'M') && (
                            (self.finish_ASAP && last_addition != '=' && !could_prevent_cycle_completion) ||
                            (self.openCycles.len() > 0 && could_prevent_cycle_completion ) ||
                            (last_addition == '=' && self.nestingCycleToClose[self.nestingCycleToClose.len() -1 ] != 0) ||
                            (last_addition == '=' && self.open_nesting_ASAP)
                        )){




                        //liaisons interdites
                        let mut prev_atom = self.nestingOpenAtom[self.nestingOpenCovalence.len()-1];
                        let mut n = 0;
                        while prev_atom == 'X' {
                            n += 1;
                            prev_atom = self.nestingOpenAtom[self.nestingOpenCovalence.len()-1 - n];
                        }

                        let mut bondType  = 1;
                        if last_addition == '=' || (last_addition == '=' && last_last_addition == '('){
                            bondType = 2;
                        }

                        if i == 'O' {
                            //println!("quoi ??? {} {}", prev_atom, bondType)
                        }

                        if legalBonds.contains(&(i, prev_atom, bondType)) || legalBonds.contains(&(prev_atom, i, bondType) )|| i == 'U' || i == 'M' || i == 'L' || i == 'W'{
                            //println!("coup legal");
                            let mv = Move{atom : i, doubleLink : false, nesting : false, closeNesting : false, cycle : 0 };
                            vec.push(mv);
                        }
                    }
                }
            }
            /*
            //ouverture d'un nesting, interdit en finish ASAP
            if  last_addition != '(' && !NUMBERS.contains(&last_addition) &&  !self.finish_ASAP {
                let mv = Move{atom : ' ', doubleLink : false, nesting : true, closeNesting : false, cycle : 0 };
                vec.push(mv);
            }
            if self.open_nesting_ASAP && ATOMS.contains(&last_addition)  {
                //rien, c'est foireux ça prend pas les nesting en compte
                let mv = Move{atom : ' ', doubleLink : false, nesting : true, closeNesting : false, cycle : 0 };

                return vec![mv];
            }

             */


        }

        // double lien et ouverture de cycle interdit en finish ASAP
        if (self.nestingOpenCovalence[self.nestingOpenCovalence.len() -1] >= 2 || (last_addition == '(' && self.nestingOpenCovalence[self.nestingOpenCovalence.len() -2] >= 1) ) && !self.finish_ASAP {
            //double lien, interdit après un départ de cycle
            if last_addition != '=' {
                let mv = Move { atom: ' ', doubleLink : true, nesting: false, closeNesting: false, cycle: 0 };
                vec.push(mv);
            }



            //ouverture d'un cycle, max 9
            // && last_addition != '='
            if self.openCycles.len() + self.closedCycles < 9  && last_addition != '(' && !NUMBERS.contains(&last_addition) {//à la base ça s'arrêtait au '='
                let mv = Move{atom : ' ', doubleLink : false, nesting : false, closeNesting : false, cycle : self.openCycles.len() + self.closedCycles + 1 };
                vec.push(mv);
            }

        }

        //println!("{:?} et heuuu {}", &self.SMILE, vec.len());

        if PRUNING_MOVE_TR >= 0.0 {
            let vec2 = vec.clone();
            let vecbackup = vec.clone();
            vec.clear();
            //println!("et bah {}", vec.len());
            for m in vec2 {
                let heuri = self.heuristic(m);
                if  heuri != 0.0 && heuri.exp() > PRUNING_MOVE_TR {
                    vec.push(m);
                }
            }


            if vec.len() == 0 && vecbackup.len() > 0 {
                //println!("NGRAM_LEN too high, no corresponding stats");
                vec = vecbackup;
            }


        }


        //println!("et bah {}", vec.len());


        return vec;
    }


    pub fn terminal(&mut self) -> bool{
        //return self.SMILE.len() > 10;

        return self.legal_moves().len() == 0;
    }

    pub fn score(&mut self) -> f64 {

        //println!("{:?}", self.SMILE);

        /*
        let mut s = String::from("");
        for c in &self.SMILE {
            s.push(*c);
        }
        println!("{}", s);

         */
        //let sc =Self::lipinskiness(self.SMILE.clone());

        if self.openCycles.len() != 0 {
            //println!("nonono {:?}", self.SMILE); //should not happen often
            //let a = vec![0]; a[2]; //debug
            return 0.0}
        if self.nestingOpenCovalence.len() != 1 {
            println!("ninini {:?}", self.SMILE); //should never happen
            return 0.0}

        let sc =self.lipinskiness(self.SMILE.clone());

        if sc >= Self::BEST_POSSIBLE_SCORE {
            self.reached_best_score = true;
        }

        return  sc; //Self::ratioH(self.SMILE.clone());  //self.SMILE.len() as f64;

    }

    pub fn backtrackCycle(SMILE : Vec<char>, looking_for : char) -> (usize, bool, Vec<usize>, usize, usize){ //return longueur, aromatique, nombre de liaisons simples devenues rigides, nombre d'atomes propres (pas contenus dans une autre boucle)
        let mut cycle_length = 0;
        let mut left_current_nesting_level = 0;//ne compte les atomes que si égal à 0
        let mut right_current_nesting_level = 0;

        let mut left_chain = Vec::new();
        let mut right_chain = Vec::new();
        let mut left_chain_id = Vec::new();
        let mut right_chain_id = Vec::new();
        let mut would_be_made_rigid = Vec::new(); //indices des atome desquels partiraient une liaison devenue rigide par aromaticité

        let mut cycle_encountered = false;

        let mut inproper_atoms = 0;
        let mut active_subcycles = vec![];

        let mut bond_cost = 1;


        let mut s = String::from("");
        for c in &SMILE {
            s.push(*c);
        }
        //println!("{} backtrack with {}", s, looking_for);


        for i in 0..SMILE.len() {
            let mut indice = SMILE.len()-1 - i;

            if SMILE[indice] == looking_for {
                cycle_encountered = true;

                if SMILE[indice-1] == '='{
                    bond_cost = 2;
                }
            }

            if NUMBERS.contains(&SMILE[indice]){
                if active_subcycles.contains(&SMILE[indice]){
                    active_subcycles.pop();
                }else{
                    active_subcycles.push(*&SMILE[indice]);
                }
            }

            if SMILE[indice] == ')'{
                right_current_nesting_level -= 1;
                if cycle_encountered {
                    left_current_nesting_level -= 1;
                }
            }

            if SMILE[indice] == '('{
                right_current_nesting_level += 1;
                if right_current_nesting_level > 0 {right_current_nesting_level = 0;}
                if cycle_encountered {
                    left_current_nesting_level += 1;
                    if left_current_nesting_level > 0 {left_current_nesting_level = 0;}
                }
            }

            if ATOMS.contains(&SMILE[indice]) {

                if left_current_nesting_level == 0 && right_current_nesting_level == 0 && cycle_encountered {
                    //fin



                    //println!("{}", cycle_length + 1);

                    let mut chain = right_chain.clone();
                    let mut chain_id = right_chain_id.clone();
                    for i in 0..left_chain.len(){
                        chain.push(left_chain[left_chain.len()-1-i]);
                        chain_id.push(left_chain_id[left_chain_id.len()-1-i]);
                    }

                    /*
                    let mut s = String::from("");
                    for c in &chain {
                        s.push(*c);
                    }
                    println!("chaine {}", s);

                     */




                    let mut aromatic = true;
                    //pas deux = directement après et nombre de = à la moitié de la chaine
                    let mut last_equal =1; //probleme : avec 1, ç a ne prend en compte que ceux de la forme C=CC=CC=, doit prendre en compte les forme =CC=CC=C
                    if chain.len() != 0 && chain[0] == '=' {last_equal = 0} //ça suffira ?

                    let mut double_count = 0;
                    let mut atom_count = 1; //on ajouts l'atome initial qui n'est aps compté
                    for i in 0..chain.len(){
                        let c = chain[i];
                        if c == '=' {
                            would_be_made_rigid.pop();

                            double_count += 1;
                            if last_equal < 2 {
                                aromatic = false;
                            }

                            if last_equal > 2 {
                                aromatic = false;
                            }

                            last_equal = 0;
                        }
                        if ATOMS.contains(&c) {
                            last_equal += 1;
                            atom_count += 1;
                            would_be_made_rigid.push(chain_id[i]);
                        }
                    }

                    if double_count*2 < atom_count-1 {
                        aromatic = false
                    }

                    let mut aromatic_rigid = atom_count - double_count;
                    if !aromatic {
                        aromatic_rigid = 0;
                        would_be_made_rigid = Vec::new();
                    }


                    return (cycle_length + 1, aromatic, would_be_made_rigid, cycle_length - inproper_atoms, bond_cost) //avant c'était aromatic rigid
                }

                if left_current_nesting_level == 0 && cycle_encountered {
                    left_chain.push(SMILE[indice]);
                    left_chain_id.push(indice);
                    cycle_length += 1;
                    if active_subcycles.len() > 0{
                        inproper_atoms += 1;
                    }
                }

                if right_current_nesting_level == 0 {
                    right_chain.push(SMILE[indice]);
                    right_chain_id.push(indice);
                    cycle_length += 1;
                    if active_subcycles.len() > 0{
                        inproper_atoms += 1;
                    }
                }
            }

            if SMILE[indice] == '='{


                if left_current_nesting_level == 0 && cycle_encountered {
                    left_chain.push(SMILE[indice]);
                    left_chain_id.push(indice);
                }

                if right_current_nesting_level == 0 {
                    right_chain.push(SMILE[indice]);
                    right_chain_id.push(indice);
                }
            }
        }



        //on a un problème si ça arrive
        println!("error, backtrack cycle should not return here");
        would_be_made_rigid = Vec::new();
        return (cycle_length, false, would_be_made_rigid, cycle_length - inproper_atoms, bond_cost);
    }

    pub fn ratioH(SMILE : Vec<char>) -> f64{
        let mut hydrog = 2.0;
        let mut other_atoms = 0.0;
        for c in SMILE {
            if c == 'C' {hydrog += 2.0; other_atoms += 1.0;} //liasons -2
            if c == 'O' {hydrog += 0.0; other_atoms += 1.0;}
            if c == 'N' {hydrog += 1.0; other_atoms += 1.0;}
            if c == 'F' {hydrog += -1.0; other_atoms += 1.0;}
            if c == 'S' {hydrog += 4.0; other_atoms += 1.0;}
            if c == 'U' {hydrog += 0.0; other_atoms += 3.0;}
            if c == 'M' {hydrog += 0.0; other_atoms += 2.0;}
            if c == 'L' {hydrog += -1.0; other_atoms += 1.0;}
            if c == 'W' {hydrog += -1.0; other_atoms += 1.0;}
            if c == 'P' {hydrog += 3.0; other_atoms += 1.0;}
            if c == '=' {hydrog += -1.0}
            if NUMBERS.contains(&c) {hydrog += -1.0}

        }
        return other_atoms/hydrog;//hydrog/other_atoms;
        //return hydrog/other_atoms;
    }

    pub fn molecularWeight(SMILE : Vec<char>)-> f64{
        let mut mass = 2.0;
        //S(=O)(=O) -> U
//S(=O) -> M
//Cl -> L
//C(F)(F)(F) -> W
        for c in SMILE {
            if c == 'C' {mass += 14.0} //masse atomique + liasons -2
            if c == 'O' {mass += 16.0}
            if c == 'N' {mass += 15.0}
            if c == 'F' {mass += 18.0}
            if c == 'S' {mass += 32.0}
            if c == 'P' {mass += 32.0}

            if c == 'U' {mass += 64.0}
            if c == 'M' {mass += 48.0}
            if c == 'L' {mass += 34.0}
            if c == 'W' {mass += 68.0}

            if c == '=' {mass += -1.0}
            if NUMBERS.contains(&c) {mass += -1.0}

        }
        return mass;
    }

    pub fn make_from_string(s : &str) -> State{

        let mut st = State::new();
        for c in s.chars() {
            // do something with `c`

            let mut m : Move = Move { atom: ' ', doubleLink : false, nesting: false, closeNesting: false, cycle: 0 };

            if c == '(' {m.nesting = true;}
            if c == ')' {m.closeNesting = true;}
            if c == '=' {m.doubleLink = true;}
            if NUMBERS.contains(&c) {m.cycle = c.to_digit(10).unwrap() as usize;}
            if ATOMS.contains(&c) {m.atom = c;}
            st.play(m);
        }
        return st;
    }

    pub fn lipinskiness(&self, SMILE : Vec<char>) -> f64{
        let mut total = 0.0;

        //nitrogen, oxygen, fluor, ne peut donner que s'il possède un hydrogène, mais peut toujours accepter
        let mut n_hbond_donor = 0.0;
        let mut n_hbond_acceptor = 0.0;
        for c in 0..SMILE.len() {
            let mut possede_hydro = 0;
            if c > 0 && SMILE[c-1] == '=' {
                possede_hydro -= 1
            }
            if SMILE[c] == 'O' {possede_hydro += 1}
            if SMILE[c] == 'N' {possede_hydro += 2}

            let mut current_level = 0;
            for c2 in c+1..SMILE.len(){
                if possede_hydro == 0 {
                    break
                }
                if SMILE[c2] == ')' {
                    current_level -=1;
                }
                if SMILE[c2] == ')' {
                    current_level +=1;
                }
                if current_level == 0{
                    if SMILE[c2] != ')' && SMILE[c2] != '(' {
                        possede_hydro -=1;
                    }
                }
            }

            if possede_hydro > 0 {
                if SMILE[c] == 'O' {n_hbond_donor += 1.0}
                if SMILE[c] == 'N' {n_hbond_donor += 1.0}
            }else{
                if SMILE[c] == 'O' {n_hbond_acceptor += 1.0}
            }

            if SMILE[c] == 'N' {n_hbond_acceptor += 1.0}
            if SMILE[c] == 'F' {n_hbond_acceptor += 1.0}//F est par nature accepteur uniquement
            if SMILE[c] == 'W' {n_hbond_acceptor += 3.0}//F est par nature accepteur uniquement
            if SMILE[c] == 'U' {n_hbond_acceptor += 2.0}
            if SMILE[c] == 'M' {n_hbond_acceptor += 1.0}
        }

        let mol_w = Self::molecularWeight(SMILE.clone());

        let mut n_rings = 0.0;
        for &c in &SMILE {
            if NUMBERS.contains(&c) {n_rings += 0.5}
        }


        let mut hydrog = 2.0;
        let mut heavy_atoms = 0.0;
        let mut carbon_count = 0.0;
        let mut nitro_count = 0.0;
        let mut oxy_count = 0.0;
        let mut sulfur_count = 0.0;
        for &c in &SMILE {
            if c == 'C' {hydrog += 2.0; heavy_atoms += 1.0; carbon_count += 1.0;} //liasons -2
            if c == 'O' {hydrog += 0.0; heavy_atoms += 1.0; oxy_count += 1.0;}
            if c == 'N' {hydrog += 1.0; heavy_atoms += 1.0; nitro_count += 1.0;}
            if c == 'F' {hydrog += -1.0; heavy_atoms += 1.0;}
            if c == 'S' {hydrog += 0.0; heavy_atoms += 1.0; sulfur_count += 1.0;}
            if c == 'P' {hydrog += 3.0; heavy_atoms += 1.0;}

            if c == 'U' {hydrog += 0.0; heavy_atoms += 3.0; sulfur_count += 1.0; ; oxy_count += 2.0;}
            if c == 'L' {hydrog += -1.0; heavy_atoms += 1.0;}
            if c == 'M' {hydrog += 0.0; heavy_atoms += 2.0; sulfur_count += 1.0; ; oxy_count += 1.0;}
            if c == 'W' {hydrog += -1.0; heavy_atoms += 4.0;}

            if c == '=' {hydrog += -1.0}
            if NUMBERS.contains(&c) {hydrog += -1.0}

        }
        //println!("hydrog : {} others : {}", hydrog, other_atoms);
        let mut n_atoms = hydrog + heavy_atoms;




        let mut n_rigid_bonds = 0.0; //double liaisons et liaisons dans cycles aromatiques


        let mut SMILE_sofar = Vec::new();
        let mut alreadyEncountered = Vec::new();
        let mut rigid_bonds = HashSet::new();

        let mut ring_sizes = Vec::new();
        let mut aromatic_rings = 0.0;
        for &c in &SMILE{
            if c == '=' {
                n_rigid_bonds += 1.0;
            }
            if NUMBERS.contains(&c) {
                if alreadyEncountered.contains(&c) {
                    let (s, a, ri, pro, bc) = Self::backtrackCycle(SMILE_sofar.clone(), c); //TODO : vérifier que ça calcule comme il faut
                    ring_sizes.push(s);
                    if a {aromatic_rings +=1.0};
                    //println!("cycle {} donne {} liaisons rigides", c, ri.len());
                    //println!("rigidified : {:?}", ri);
                    //n_rigid_bonds += ri.len() as f64;
                    for i in ri {
                        rigid_bonds.insert(i);
                    }

                }else{
                    alreadyEncountered.push(c);
                }
            }
            SMILE_sofar.push(c);
        }
        n_rigid_bonds += rigid_bonds.len() as f64;

        let mut n_rotatable_bonds = 0.0; //les autres, taille de la molécule - 1 - n_rigir_bonds
        let mut all_bonds = 0.0;

        for &c in &SMILE {
            all_bonds += 1.0;
            if NUMBERS.contains(&c) || c == '(' || c == ')' {
                all_bonds += -0.5;
            }
        }
        n_rotatable_bonds = all_bonds - n_rigid_bonds;


        let mut n_ring_size_six = 0.0;
        for i in ring_sizes { if i == 6 {n_ring_size_six += 1.0}}


        //score += -std::cmp::max(n_hbond_donor-5.0,0 )/5.0;

        let mut all_attributes = Vec::new();
        all_attributes.push(-Self::max(n_hbond_donor-5.0,0.0 )/5.0);
        all_attributes.push(-Self::max(n_hbond_acceptor-10.0,0.0 )/10.0);

        //all_attributes.push(Self::min(n_rings,3.0 )/3.0); //au moins trois anneaux
        //all_attributes.push(1.0);
        //all_attributes.push(Self::min(n_rigid_bonds,18.0 )/18.0); //absolu pas bien
        //all_attributes.push(Self::min(n_rotatable_bonds,6.0 )/6.0);

        //grand
        //all_attributes.push(Self::min(n_atoms,50.0 )/50.0);
        //all_attributes.push(-Self::max(n_atoms-100.0,0.0 )/100.0);
        //all_attributes.push(-Self::max(mol_w-1000.0,0.0 )/1000.0);

        //moyen (lipinski)
        all_attributes.push(Self::min(n_atoms,20.0 )/20.0);
        //all_attributes.push(-Self::max(20.0-n_atoms,0.0 )/20.0); //tentative d'éviter le minimum local pour UCT (bof)
        all_attributes.push(-Self::max(n_atoms-70.0,0.0 )/70.0);
        all_attributes.push(-Self::max(mol_w-500.0,0.0 )/500.0);

        //plus petit
        //all_attributes.push(Self::min(n_atoms,20.0 )/20.0);
        //all_attributes.push(-Self::max(n_atoms-30.0,0.0 )/30.0);
        //all_attributes.push(-Self::max(mol_w-200.0,0.0 )/200.0);

        all_attributes.push(Self::min(n_rigid_bonds,n_atoms/6.0 )/(n_atoms/6.0)); //au moins 1 rigid bond pour 6 atomes
        //all_attributes.push(-Self::max(n_rigid_bonds-n_atoms/2.0,0.0 )/(n_atoms/2.0));//au plus 3 rigid bond pour 6 atomes

        //mes ajouts
        all_attributes.push(-Self::max(n_rings - 5.0,0.0 )/5.0); //max 5 anneaux stp
        //all_attributes.push(Self::min(n_ring_size_six,2.0 )/2.0); //au moins deux benzene

        //all_attributes.push(-Self::max(sulfur_count*10.0-carbon_count,0.0 )/carbon_count); //au maximum un sulfur pour 10 carbones
        //all_attributes.push(-Self::max(oxy_count*2.0-carbon_count,0.0 )/carbon_count); //au maximum un oxy pour 2 carbones
        //all_attributes.push(-Self::max(nitro_count*2.0-carbon_count,0.0 )/carbon_count); //au maximum un nitro pour 2 carbones


        //if n_rings >= 2.0 {all_attributes.push(-Self::max(n_rings-1.0-aromatic_rings,0.0 ));} //si au moins deux cycle alors aromatique (cause à ne faire qu'un cycle)
        //all_attributes.push(-Self::max(n_rings-aromatic_rings,0.0 )); //tous anneaux aromatiques (cause à ne pas faire de cycle)


        if self.target_NtoC_ratio >= 0.0 {
            let NtoC_ratio = nitro_count/carbon_count;
            all_attributes.push(-Self::max((self.target_NtoC_ratio-NtoC_ratio).abs() - 0.1,0.0 ));
        }

        if self.target_OtoC_ratio >= 0.0 {
            let OtoC_ratio = oxy_count/carbon_count;
            all_attributes.push(-Self::max((self.target_OtoC_ratio-OtoC_ratio).abs() - 0.1,0.0 ));
        }


        for i in &all_attributes {
            total += i;
        }

        //println!(" attributs: {:?}", all_attributes);

        return 1000.0 + total;

    }

    pub fn min(a : f64, b : f64) -> f64{
        if a < b {
            return a;
        }
        return b;
    }
    pub fn max(a : f64, b : f64) -> f64{
        if a > b {
            return a;
        }
        return b;
    }

    pub fn heuristic(&mut self, m : Move) -> f64{

        //println!("appel heuri");

        let mut mv = 'a';
        if m.doubleLink {mv = '='}
        if m.nesting {mv = '('}
        if m.closeNesting {mv = ')'}
        if m.cycle != 0 {mv = 'X'}
        if m.atom != ' ' {mv = m.atom}

        if HEURISTIC_MODE == "neural" {

            let mut SMILEstring: Vec<String> = Vec::new();
            for c in &self.SMILE {

                SMILEstring.push(c.to_string())
            }

            let mut SMILEstr: Vec<&str> = Vec::new();
            for c in &SMILEstring {

                SMILEstr.push(&c.as_str())
            }




            if self.storedPrior.label.len() == 0 {
                //println!("costly"); // 0.056s par execution du réseau
                self.storedPrior = prior.predict(SMILEstr).unwrap();
            }else{
                //println!("avoided a prediction");
            }
            let mut pred = self.storedPrior.clone();


            //exact same NN as chemTS
            //let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "c", "1", "2", "o", "=", "O", "N", "3", "F", "[C@@H]", "n", "-", "#", "S", "Cl", "[O-]", "[C@H]", "[NH+]", "[C@]", "s", "Br", "/", "[nH]", "[NH3+]", "4", "[NH2+]", "[C@@]", "[N+]", "[nH+]", "\\", "[S@]", "5", "[N-]", "[n+]", "[S@@]", "[S-]", "6", "7", "I", "[n-]", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[o+]", "[CH2-]", "[CH-]", "[SH+]", "[O+]", "[s+]", "[PH+]", "[PH]", "8", "[S@@+]"];

            //same but with explicit cycles
            //let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", "S", "Cl", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "4", "[NH2+]", "[C@@]", "[N+]", "\\", "[S@]", "5", "[N-]", "[S@@]", "[S-]", "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", "[S@@+]"];
            //same but with explicit cycles and shortcuts
            let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", "S", "L", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "W", "4", "[NH2+]", "U", "[C@@]", "[N+]", "\\", "M", "[S@]", "5", "[N-]", "[S@@]", "[S-]", "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", "[S@@+]"];

            //return pred.unwrap().confidence as f64;

            //in case of loop
            if mv == 'X'{
                let mut sum_conf = 0.0;
                for i in 0..pred.label.len(){
                    if val[pred.label[i] as usize].parse::<f64>().is_ok() {
                        sum_conf += pred.confidence[i] as f64;
                    }
                }
                return sum_conf.ln();
            }

            //in case of atom
            if mv == m.atom {
                let mut sum_conf = 0.0;
                for i in 0..pred.label.len(){
                    if val[pred.label[i] as usize] == mv.to_string().as_str() || val[pred.label[i] as usize] == mv.to_lowercase().to_string().as_str()  {
                        sum_conf += pred.confidence[i] as f64;
                    }
                }
                return sum_conf.ln();
            }

            //in case of anything else
            for i in 0..pred.label.len(){
                if val[pred.label[i] as usize] == mv.to_string().as_str()  {
                    return pred.confidence[i].ln() as f64;
                }
            }

        }

        if HEURISTIC_MODE == "ngram" {
            let mut s = String::from("");
            //let mut canCloseCycle = false;
            //prend en compte la longueur du cycle
            if mv == 'X' && self.openCycles.contains(&m.cycle){
                //canCloseCycle = true;
                let (si, a, ri, pro, bc) = Self::backtrackCycle(self.SMILE.clone(), char::from_digit(m.cycle as u32, 10).unwrap());

                /*
                nombre de taille 5 : 0.595046
                nombre de taille 6 : 1.958514
                nombre de taille 7 : 0.053870
                 */
                if si == 5 { return (0.595046f64/(0.595046f64+1.958514f64+0.053870f64)).ln(); }
                if si == 6 { return (1.958514f64/(1.958514f64+0.053870f64)).ln(); }
                if si == 7 { return (1.0f64).ln(); }

                /*
                if s < 10 {mv = char::from_digit(s as u32, 10).unwrap()}else{
                    mv = '0';
                }
                */
            }

            for i in 0..NGRAM_LEN {
                if(self.SMILE.len() > NGRAM_LEN -i-1){
                    if NUMBERS.contains(&self.SMILE[self.SMILE.len() - (NGRAM_LEN -i)]) {
                        s.push('X');
                    }else{
                        s.push(self.SMILE[self.SMILE.len()-(NGRAM_LEN -i)]);
                    }

                }
            }


            /*

            let mut sum_close_cycle = 0.0;
            for i in  0..10 {
                let key = ngrams.get(&*s);
                if key.is_none() {
                    //skip
                }else{
                    let ret = key.unwrap().nextAtom.get(&char::from_digit(i as u32, 10).unwrap());
                    if ret.is_none() {
                        //skip
                    }else{
                        sum_close_cycle += ret.unwrap();
                    }
                }
            }
            */




            let key = ngrams.get(&*s);
            if key.is_none() {
                //println!("waaaa {}", s);
                return 0.0
            }
            //println!("{}", s);
            //println!("hey {:?}", key.unwrap().nextAtom);


            let ret = key.unwrap().nextAtom.get(&mv);
            if ret.is_none() {
                return 0.0
            }



            /*
            if sum_close_cycle != 0.0 && sum_close_cycle != 1.0 {
                if canCloseCycle {return (*ret.unwrap()/sum_close_cycle).ln();}
                return (*ret.unwrap()*(1.0-sum_close_cycle)).ln();
            }

             */
            return (*ret.unwrap()).ln();



        }

        return 0.0;

    }



    pub fn smoothedScore(&mut self) ->f64{
        return self.score();
    }

}