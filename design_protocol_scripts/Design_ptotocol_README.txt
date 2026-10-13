The computational design method for designing hetero-chiral binder_delta_sap

======================== Step 1: Initialization ========================
Make a directory for this project and cd to it. Then make the following directories:
PROJECT/1.input
PROJECT/2.rifgen
PROJECT/3.scoffold_polyVal
PROJECT/4.patchdock
PROJECT/5.rifdock
PROJECT/6.pre_handle
PROJECT/7.predictor_1
PROJECT/8.pilot_1
PROJECT/9.fastdesign_1
PROJECT/10.motif_graft
PROJECT/11.predictor_2
PROJECT/12.pilot_2
PROJECT/13.fastdesign_2
PROJECT/14.final_filter

Prepare your scaffolds list.


======================== Step 2: Target Preparation ========================
Check the input folder for preparation of D-version target. You can use this(you may need the pyrosetta environment for this script):
    python invert_chiral.py infile.pdb outfile.pdb
    mv outfile.pdb target.pdb

Make target 1 chain and change it to chain A. The following should work:
    cat target.pdb | grep '^ATOM' | sed 's/./A/22' > tmp.pdb  
    mv tmp.pdb target.pdb


======================== Step 3: Rifgen ========================
You first have to select your rifgen residues, and put them in a text file called rifgen_residues.list. One at a time on seperate lines.

You will need a rifgen.flag file. The only lines that should need editing are the ones in the I/O section.

Once you are happy with your settings, you can schedule your rifgen run. The rifgen is available in Rifdock(https://github.com/rifdock/rifdock)  It is critical that you save the log. These runs often take about 8 hours. 16 cpu and 100g of ram is recommended.

After your rifgen has finished, look in the output folder for your target pdb. This is the pdb that you will use for patchdock. It is also important convert your output pdb to chainB.


======================== Step 4: Scoffold_polyVal ========================
You need to convert your scaffolds to poly-valine. You may need this python script: mutate_polyXXX.py


======================== Step 5: Patchdock ========================
Once you have your poly-valine scaffolds, it's time to setup the patchdock commands. The following may works:
    patch_dock.Linux ${scaf}.params ${scaf}.out


======================== Step 6: Rifdock ========================
To begin, you will need a rifdock flags file. (rifdock.flags) 

Next, locate your rifgen.log file. Open this file and scroll to the bottom. You will see a large sections denoted with:  what you need for docking  Copy and paste this section to the top of your rifdock.flag. 

Then you can run rifdock. The Rifdock(https://github.com/rifdock/rifdock) can be download and compiled.


======================== Step 7: Pre_handle ========================
After you get your rifdock outputs, you may change the residue name to D-version and renumber the residue index in each chain. change_L_name_to_D_name.sh can help change name, and renum_chain.pl can renumber index.


======================== Step 8: Predictor_1 ========================
RifDock is very good at finding good docks, but it is not the best at scoring docks. For this reason, a lot of your RifDock outputs are not worth designing. In this step, we'll filter out the bad rifdock outputs so you don't waste cpu time designing them.  

For this filtering, rosetta script is performed, you will need predictor.xml and predictor.flags to run the predictor with rosetta_scripts.


======================== Step 9: Pilot_1 ========================
We also need to setup some pilot runs for machine learning that we will be doing shortly. 5000 pdbs will be randomly selected.

You will need design.xml and pilot.flags to run the pilot with rosetta_scripts.


======================== Step 10: Fastdesign_1 ========================
After above two jobs have finished, it's time to collect the results. predictor_runs.sc now contains the predictor scores for all of your rifdock outputs. pilot_runs.sc now contains the predictor scores for all of choosed(5000 pdbs) rifdock outputs.

We now need to use the pilot runs to select which of your designs should move forward for FastDesign. Predictor_notebook is need in the notebook folder. Then you can running the Predictor_notebook with jupyter notebook and choose pdbs you will do the fastdesign runs.(you may choose 10000 pdbs for fastdesign_1)

For fastdesign_1, you will need design.xml and fastdesign.flags to run the fastdesign with rosetta_scripts. 


======================== Step 11: Motif_graft ========================
In this step, we will scan through your fastdesign outputs and extract secondary structure elements. The elements will be your motifs. you will need this scripts to extract your motifs: extract_PPI_motifs.py.

After those commands finish, it's time to cluster your motifs. These commands will cluster them:
    motif_clustering/cluster all_cluster.list 0.5 a (yes really the 'a' needs to be there)
    python convert_cluster_output.py cluster_results.list > cluster_results.dat

Now that your motifs are clustered, it's time to use another jupyter notebook. The notebook is here: notebook/look_at_motifs_by_length_w_hbond.ipynb. Run this notebook and it will output motifs for both helices and strand pairs.

We need to get the hotspots from your motifs. This uses anything with a ddG better than -2 as a hotspot. (When compared to alanine)
    python get_path_and_hotspot.py <(cat motifs_by_size/*.list) > motifs_to_use_w_hotspots.list

Then you will need graft.xml and graft.flags to run the grafting with rosetta_scripts.


======================== Step 12: Predictor_2 ========================
Now we have grafting output pdbs, however, there may be more than 10M outputs, so we have to filter all the outputs and select pdbs we want to design. 

The process is very similar to Predictor_1, the only different between predictor.xml and predictor2.xml is line 52(predictor2.xml). Since we didn't renumber all the output from grafting.

You will need predictor2.xml and predictor2.flags to run the predictor with rosetta_scripts.


======================== Step 13: Pilot_2 ========================
We also need to setup some pilot runs for machine learning that we will be doing shortly. 5000 pdbs will be randomly selected.

You will need design2.xml and pilot2.flags to run the pilot with rosetta_scripts.


======================== Step 14: Fastdesign_2 ========================
After above two jobs have finished, it's time to collect the results. predictor_runs.sc now contains the predictor scores for all of your rifdock outputs. pilot_runs.sc now contains the predictor scores for all of choosed(5000 pdbs) rifdock outputs.

We now need to use the pilot runs to select which of your designs should move forward for FastDesign. Predictor_notebook is need in the notebook folder. Then you can running the Predictor_notebook with jupyter notebook and choose pdbs you will do the fastdesign runs.(you may choose 100000 pdbs for fastdesign_1)

For fastdesign_2, you will need design2.xml and fastdesign2.flags to run the fastdesign with rosetta_scripts.


======================== Step 15: Final filter ========================
Now you will get youe 100000 pdbs from fastdesign_2, we will do a final filter to select pdbs with good matrics.

For final_filter, you will need filter.xml and filter.flags to run the fastdesign with rosetta_scripts. 

After above jobs finished. We will collect the results. final.sc contains all the scores of your final_filtering outputs. We will run a final_filter notebook to select pdbs we want. The notebook you need is here: notebook/final_filtering.ipynb. 

Finally, approximately 10000 pdbs will survived from the notebook and further experiment will be carried. 


