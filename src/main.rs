#![allow(non_snake_case)]

#[macro_use]
extern crate lazy_static;

use std::collections::HashMap;
use std::io;
use std::io::Write;
use std::time::Instant;
use calc::mean;
use calc::std;
use crate::tools::calc;
use crate::tools::resultSaver;

use chrono::{DateTime, Utc};

mod tools;
mod methods;
mod models;

fn main() {
    println!("Hello, world!");

    let now: DateTime<Utc> = Utc::now();
    println!("UTC now is: {}", now);


    testMolGen();





}


pub fn testMolGen(){

    let mut molGenState = models::SMILESgen::State::new();


    if false{

        //C1(C=CC=C(C pour tester les nouveaux SMILES
        molGenState = models::SMILESgen::State::make_from_string("C1=C(CCC");

        println!("open cov {:?}", molGenState.nestingOpenCovalence);
        println!("fdsfgsdf {}", molGenState.open_nesting_ASAP);
        println!("dfsfsdf {:?}", molGenState.nestingCycleToClose);
        println!("finish asap {}", molGenState.finish_ASAP);

        //let m0 = molGenState.legal_moves()[0];
        //molGenState.heuristic(m0).exp(); //initialize the NN

        //println!("label : {:?}", molGenState.storedPrior.label);
        //println!("prior : {:?}", molGenState.storedPrior.confidence);

        let start_time: Instant = Instant::now();

        for m in molGenState.legal_moves() {
            println!("atom : {}, nesting : {}, close nesting : {}, cycle : {}, double : {}  ", m.atom, m.nesting, m.closeNesting, m.cycle, m.doubleLink);
            println!("prior value {}", molGenState.heuristic(m).exp());
        }
        println!("len : {}", molGenState.legal_moves().len());
        println!("time : {} ",  start_time.elapsed().as_secs_f64());

        while true {}
    }

    let mut v = true;

    let FDA_OtoC_ratio = vec![0.117, 0.135, 0.156, 0.139, 0.113, 0.115, 0.061, 0.037, 0.045, 0.009, 0.016, 0.013, 0.005, 0.006, 0.002, 0.003, 0.003, 0.001, 0.001, 0.013];
    let FDA_NtoC_ratio = vec![0.266, 0.180, 0.184, 0.128, 0.076, 0.058, 0.031, 0.011, 0.016, 0.003, 0.016, 0.011, 0.001, 0.003, 0.000, 0.001, 0.002, 0.000, 0.000, 0.004];
    let mut OtoC_bins = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut NtoC_bins = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut OtoC_ratio = OtoC_bins.clone();
    let mut NtoC_ratio = NtoC_bins.clone();

    let mut molecule_count = 0;

    let mut start_time = Instant::now();

    while molecule_count < 1000{

        let mut OtoC_diff = 0.0;
        let mut target_OtoC_ratio = 0.0;

        for i in 0..FDA_OtoC_ratio.len(){
            if FDA_OtoC_ratio[i] - OtoC_ratio[i] > OtoC_diff {
                OtoC_diff = FDA_OtoC_ratio[i] - OtoC_ratio[i];
                target_OtoC_ratio = (i as f64)/(FDA_OtoC_ratio.len() as f64)
            }
        }

        let mut NtoC_diff = 0.0;
        let mut target_NtoC_ratio = 0.0;

        for i in 0..FDA_NtoC_ratio.len(){
            if FDA_NtoC_ratio[i] - NtoC_ratio[i] > NtoC_diff {
                NtoC_diff = FDA_NtoC_ratio[i] - NtoC_ratio[i];
                target_NtoC_ratio = (i as f64)/(FDA_NtoC_ratio.len() as f64)
            }
        }


        let mut targetState = molGenState.clone();

        //uncomment here to enable ratio enforcement
        //targetState.target_NtoC_ratio = target_NtoC_ratio;
        //targetState.target_OtoC_ratio = target_OtoC_ratio;





        let mut st = methods::NMCS::launch_nmcs(targetState,3, 1.0, v,0.0, String::from("NMCS_SMILEGEN"));
        //let mut st = methods::Sampling::Sample(targetState,  1000000000, 1.0, 0.0, String::from("Sampling_SMILEGEN"));
        //let mut st = methods::UCT::launch_UCT(targetState, 1.0, 1000000000, 0.0, 10.0, String::from("UCT_SMILEGEN"));
        //let mut st = methods::PUCT::launch_PUCT(targetState, 1.0, 1000000000, 1.0, 0.0, String::from("PUCT_SMILEGEN"));

        if st.reached_best_score {
            molecule_count += 1;
            //let st = methods::NMCS::launch_nmcs(molGenState.clone(),3, 0.0, v,0.0, String::from("NMCS_network"));
            //println!("meilleur score : {}", st.score());
            //println!("{:?}", st.SMILE);

            let mut s = vec![];
            for c in &st.SMILE {
                s.push(*c);

                if s.len() > 3 && s[s.len()-1] == ')' && s[s.len()-3] == '(' {
                    if s[s.len()-2] == '1' || s[s.len()-2] == '2' || s[s.len()-2] == '3' || s[s.len()-2] == '4' || s[s.len()-2] == '5' || s[s.len()-2] == '6' || s[s.len()-2] == '7' || s[s.len()-2] == '8' || s[s.len()-2] == '9' {
                        s.remove(s.len()-1);
                        s.remove(s.len()-2);
                    }

                }
            }

            let mut s2 = String::from("");

            for c in &s {
                if c == &'U' {
                    s2 = s2+"S(=O)(=O)";
                }

                if c == &'W' {
                    s2 = s2+"C(F)(F)(F)";
                }

                if c == &'L' {
                    s2 = s2+"Cl";
                }

                if c == &'M' {
                    s2 = s2+"S(=O)";
                }

                if c != &'U' && c != &'W' && c != &'L' && c != &'M' {
                    s2.push(*c);
                }

            }
            println!("{}", s2);
            resultSaver::writeLine(s2+"\n", String::from("SMILES generated/testsCycles/tests2"));
            //println!("{} with targets ratio O : {} N : {} ", s, target_OtoC_ratio, target_NtoC_ratio);


            let mut n_carb = 0.0;
            let mut n_nitro = 0.0;
            let mut n_oxy = 0.0;
            for &c in &st.SMILE {
                if c == 'C' {n_carb += 1.0}
                if c == 'N' {n_nitro += 1.0}
                if c == 'O' {n_oxy += 1.0}
            }

            for i in 0..OtoC_bins.len() {
                if n_oxy/n_carb < (i as f64+1.0)/OtoC_bins.len() as f64 {
                    OtoC_bins[i] += 1.0;
                    break
                }
            }

            for i in 0..NtoC_bins.len() {
                if n_nitro/n_carb < (i as f64+1.0)/NtoC_bins.len() as f64 {
                    NtoC_bins[i] += 1.0;
                    break
                }
            }

            let mut sum = 0.0;
            for i in 0..OtoC_bins.len(){
                sum += OtoC_bins[i];
            }
            for i in 0..OtoC_ratio.len(){
                OtoC_ratio[i] = OtoC_bins[i]/sum
            }

            sum = 0.0;
            for i in 0..NtoC_bins.len(){
                sum += NtoC_bins[i];
            }
            for i in 0..NtoC_ratio.len(){
                NtoC_ratio[i] = NtoC_bins[i]/sum
            }
            //println!("{:?}", NtoC_ratio);
            v =false;
        }else{
            //println!("timeout");
        }


    }

    println!("Total time : {}", start_time.elapsed().as_secs_f64());




}

