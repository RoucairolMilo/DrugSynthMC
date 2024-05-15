use anyhow::Result;
use serde::Serialize;
use tensorflow::{
    Graph, Operation, SavedModelBundle, SessionOptions, SessionRunArgs, Tensor,
    DEFAULT_SERVING_SIGNATURE_DEF_KEY, PREDICT_INPUTS, PREDICT_OUTPUTS,
};

#[derive(Debug, Serialize, Clone)]
pub struct Prediction {
    pub label: Vec<u8>,
    pub confidence: Vec<f32>,
}

pub struct Model {
    bundle: SavedModelBundle,
    input_op: Operation,
    input_index: i32,
    output_op: Operation,
    output_index: i32,
}

impl Model {
    pub fn from_dir(export_dir: &str) -> Result<Self> {
        const MODEL_TAG: &str = "serve";
        let mut graph = Graph::new();
        let bundle =
            SavedModelBundle::load(&SessionOptions::new(), &[MODEL_TAG], &mut graph, export_dir)?;
        let sig = bundle
            .meta_graph_def()
            .get_signature(DEFAULT_SERVING_SIGNATURE_DEF_KEY)?;

        let input_info = sig.get_input("embedding_input")?;
        let output_info = sig.get_output("time_distributed")?;
        let input_op = graph.operation_by_name_required(&input_info.name().name)?;
        let output_op = graph.operation_by_name_required(&output_info.name().name)?;
        let input_index = input_info.name().index;
        let output_index = output_info.name().index;

        println!("ok ?");

        Ok(Self {
            bundle,
            input_op,
            input_index,
            output_op,
            output_index,
        })
    }

    pub fn predict(&self, SMILES: Vec<&str>) -> Result<Prediction> {

        //val = [ "\n", "&","C", "(", ")", "c", "1", "2", "o", "=", "O", "N", "3", "F", "[C@@H]", "n", "-", "#", "S", "Cl", "[O-]", "[C@H]", "[NH+]", "[C@]", "s", "Br", "/", "[nH]", "[NH3+]", "4", "[NH2+]", "[C@@]", "[N+]", "[nH+]", "\\", "[S@]", "5", "[N-]", "[n+]", "[S@@]", "[S-]", "6", "7", "I", "[n-]", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[o+]", "[CH2-]", "[CH-]", "[SH+]", "[O+]", "[s+]", "[PH+]", "[PH]", "8", "[S@@+]"]

        //exact same NN as chemTS
        //let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "c", "1", "2", "o", "=", "O", "N", "3", "F", "[C@@H]", "n", "-", "#", "S", "Cl", "[O-]", "[C@H]", "[NH+]", "[C@]", "s", "Br", "/", "[nH]", "[NH3+]", "4", "[NH2+]", "[C@@]", "[N+]", "[nH+]", "\\", "[S@]", "5", "[N-]", "[n+]", "[S@@]", "[S-]", "6", "7", "I", "[n-]", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[o+]", "[CH2-]", "[CH-]", "[SH+]", "[O+]", "[s+]", "[PH+]", "[PH]", "8", "[S@@+]"];

        //same but with explicit cycles
        //let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", "S", "Cl", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "4", "[NH2+]", "[C@@]", "[N+]", "\\", "[S@]", "5", "[N-]", "[S@@]", "[S-]", "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", "[S@@+]"];

        //same but with explicit cycles and shortcuts
        //let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", "S", "Cl", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "W", "4", "[NH2+]", "U", "[C@@]", "[N+]", "\\", "M", "[S@]", "5", "[N-]", "[S@@]", "[S-]", "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", "[S@@+]"];
        let val : Vec<&str> = vec!["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", "S", "L", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "W", "4", "[NH2+]", "U", "[C@@]", "[N+]", "\\", "M", "[S@]", "5", "[N-]", "[S@@]", "[S-]", "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", "[S@@+]"];


        let mut get_int : Vec<f32> = vec![1.0];

        for c in SMILES.clone()  {
            for i in 0..val.len(){
                if c == val[i]{
                    get_int.push(i as f32);
                }
            }
        }

        while get_int.len() != 81 {
            get_int.push(0.0);
        }


        const INPUT_DIMS: &[u64] = &[1,81];


        let input_tensor = Tensor::<f32>::new(INPUT_DIMS).with_values(&get_int)?;


        let mut run_args = SessionRunArgs::new();
        run_args.add_feed(&self.input_op, self.input_index, &input_tensor);
        let output_fetch = run_args.request_fetch(&self.output_op, self.output_index);
        //println!("tensor {:?}", input_tensor);
        self.bundle.session.run(&mut run_args)?;



        let output = run_args.fetch::<f32>(output_fetch)?;

        let mut confidence = vec![];
        let mut label =  vec![];
        let mut max_lab = 0;
        let mut max_conf = 0.0;
        for i in 0..output.dims()[2] {
            confidence.push(output.get(&[0, (SMILES.len()) as u64, i]));
            if output.get(&[0, (SMILES.len()) as u64, i]) > max_conf {
                max_conf = output.get(&[0, (SMILES.len()) as u64, i]);
                max_lab = i as usize;
            }

            label.push(i as u8);

        }

        //println!("conf : {:?}", confidence);
        //println!("SMILE : {:?}, char : {} {}, confidence : {}", SMILES, val[max_lab], max_lab, max_conf);

        Ok(Prediction { label, confidence })
    }
}