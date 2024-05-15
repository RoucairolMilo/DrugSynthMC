use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::Path;

pub fn writeLine(line : String, name : String ) -> std::io::Result<()>{
    if(name.len() != 0){
        let p = "results/".to_owned() + &name + ".txt";
        if Path::new(&p).exists() {
            let mut file = OpenOptions::new()
                .write(true)
                .append(true)
                .open(p)
                .unwrap();

            file.write_all(line.as_ref());
        }else{
            let mut file = File::create(p)?;
            file.write_all(line.as_ref())?;
        }
    }
    Ok(())
}