
# kb_MaSuRCA
---

A [KBase](https://kbase.us) module generated by the [KBase SDK](https://github.com/kbase/kb_sdk).


This Module was initialized with a generated example App.  To compile and run the
example App implementation, run:

    cd kb_MaSuRCA
    make          (required after making changes to $module_name.spec)
    kb-sdk test   (will require setting test user account credentials in test_local/test.cfg)

General steps to run the MaSuRCA assemblers:

    First, create a configuration file which contains the location of the compiled assembler (the 
    executable), data and assembly parameters. 
    
    Second, run the 'masurca' script which will generate from the configuration file a shell script
    'assemble.sh', which is the main driver of the assembly.

    Finally, run the script 'assemble.sh' to assemble the data.

For more help on how to modify, register and deploy the example to KBase, see the
[KBase SDK documentation](https://github.com/kbase/kb_sdk).

