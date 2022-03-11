# UQExampleBEMHeartPosition

This is an uncertainty quantification (UQ) example computing the effect of heart position on boundary element method (BEM) ECG forward computations.  This example is similar to the work of [Swenson, etal.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3362042/) but is implemented in [UncertainSCI](https://github.com/SCIInstitute/UncertainSCI) and [SCIRun](https://github.com/SCIInstitute/SCIRun).  

## Requirements

This example requires [SCIRun](https://github.com/SCIInstitute/SCIRun),  [UncertainSCI](https://github.com/SCIInstitute/UncertainSCI),  and Python3 to run.  Additionally, this example is mainly tested with unix based systems, and may require additional modification for some systems. 

## Running the Example

Modify the `run_torso_position.py` script so the  `output_dir` and `scirun_call` variables accurately reflex the relevant system.  Next, run the script in a terminal window with the command:
```bash
python3 pyscripts/run_torso_postition.py
```
The script will generate a sample set and a script to run the samples in the `nets/torso_position_model_all.srn5` network in SCIRun.  Currently, this step is not fully automated, so there will be a prompt to execute the `tmp/UQ_torso_tmp.py` in SCIRun manually.  To do this, open SCIRun and click the run script icon, which looks like a magic wand.  The advanced toolbar will need to be enable to see it (Window -> Toolbars -> Advanced).  
After clicking the run script option, navigate to the `tmp/UQ_torso_tmp.py` and select it.  This should trigger the looped execution of the network to generate the solution files needed for UncertainSCI.  Once completed, press enter in the terminal window and resume the `pyscripts/run_torso_postition.py` script.  The script should generated a sample ECG with the median solution, some sampled solutions, and quantile ranges of the soluition distribution.  

TODO:
update script with workaround to fully automate SCIRun execution. 
