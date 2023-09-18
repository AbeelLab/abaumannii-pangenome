# A comparative study of pan-genome methods for microbial organisms: *Acinetobacter baumannii* pan-genome reveals structural variation in antimicrobial resistance-carrying plasmids

Repository for the scripts used to perform a comparative study of pan-genome methods, [manuscript published in Microbial Genomics](https://doi.org/10.1099/mgen.0.000690). If you use the code, please cite:
 
```bib
@article{baumanniipangenome2021,    
  title={A comparative study of pan-genome methods for microbial organisms: Acinetobacter baumannii pan-genome reveals structural variation in antimicrobial resistance-carrying plasmids},
  author={Urhan, Aysun and Abeel, Thomas},
  journal={Microbial Genomics},
  volume={7},
  number={11},
  year={2021},
  publisher={Microbiology Society},
  doi={https://doi.org/10.1099/mgen.0.000690}
}
```

The scripts are essentially wrap-arounds to run our comparative study, the individual methods belong to their corresponding developers. So please refer to their repositories to confirm how to cite their work.
1. [Roary (v3.13.0)](http://sanger-pathogens.github.io/Roary)
2. [Ptolemy (v1.0)](https://github.com/AbeelLab/ptolemy)
3. [PPanGGoLin (v1.0.13)](https://github.com/labgem/PPanGGOLiN)
4. [PIRATE (v1.0.3)](https://github.com/SionBayliss/PIRATE)
5. [Panaroo (v1.1.2)](https://github.com/gtonkinhill/panaroo)
6. [GOATOOLS (v0.9.9)](https://github.com/tanghaibao/goatools)

'createSA.py' and 'runGOE.py' were both written in python (3.7), and the 'comparativeStudy.sh' is a bash script that includes the commands we ran to perform our comparative study.
