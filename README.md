# LINCS_project
LINCS project with Anna Budkina and Alexey Stupnikov.

These scripts allow to conduct a comparison of similar signatures search algorithms (sssa) using user's differential expression data.

To successfully launch the process you have to have the same architecture as it is in current github folder (so I would recommend to just download all the directories and put them in the exact same folder).

Sequence of launching the scripts (as input - folder with results of differential expression analysis - tables of differential expression):
1. Launch "scripts_to_obtain_signatures_and_synonyms/SIGNATURE_OBTAINER.py" - it will create a folder with files, each containing a list of DEGs (signature), corresponding to a DE experiment.
2. Launch "scripts_to_obtain_signatures_and_synonyms/SYNONYMS_FINDER.py" - it will find all possible synonyms of a chemical (perturbagene), corresponding to a DE experiment (you have to supply to it a DE_experiment/perturbation correspondence table -- examples could be found in folder "annotation_parsed")
3. Launch any of search algorithms' scripts - they can be found in "scripts_to_run_search_algorithms" folder.
4. Launch "scripts_to_generate_plots/BARPLOT_GENERATOR.R" - to visualise the result of search algorithms work.



How to launch search algorithms' scripts:

L2S2, L1000FWD and iLINCS have identical options:

compulsory option: name of inhibitors group (e.g. CDK, MEK, PI3K, mTOR)
other options:  --project_path  ---  path to a directory with the whole project
                --out_folder  ---  NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.
                --sng_folder --- NAME of a folder with all signatures, each signature in its own file (with strict-form table of DEGs)


launch example:
python3 ~/LINCS_project_folder/scripts_to_run_search_algorithms/L2S2_ALGORITHM_RUNNER.py CDK --project_path ~/LINCS_project_folder --out_folder L2S2_results_defaults --sgn_folder signatures_defaults &> ~/LINCS_project_folder/CDK_inhibitors/L2S2.logs.txt

directory architecture for search scripts to correctly process input:

```bash
LINCS_project/
|
├── annotations_parsed/                       # Parsed annotations
|
├── scripts_to_run_search_algorithms/         # Scripts to launch algorithms
│   ├── L1000_ALGORITHM_RUNNER.py             
│   ├── L2S2_ALGORITHM_RUNNER.py             
│   └── iLINCS_ALGORITHM_RUNNER.R            
│
└── <inh_group>_inhibitors/                  # Folder of the exact type of inhibitors (e.g. <inh_group> = mTOR)
    ├── mTOR_perturbagens_synonyms.json          # Json with chemical synonyms
    ├── L1000/                                   # Folder of the exact sssa (e.g. L1000FWD)
    │   ├── L1000_results_with_empty_files/          # Folder with sssa work results 
    │   └── visualization/                           # Folder with sssa work visualizatioin 
    ├── L2S2/                                    # The same for L2S2
    ├── iLINCS/                                  # The same for iLINCS
    ├── signatures/                              # FOLDER WITH SIGNATURES LIES HERE - IN THE FOLDER OF THE EXAT TYPE OF INHIBITORS
    └── visualization/
```

On the scheme above you can also see the compulsory architecture for the sssa script to correcly process its output.




