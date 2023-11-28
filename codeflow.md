```sh
entry -> driver.cli
    driver.run
        Runtime[].do_workflow

Runtime[]
    Read the yaml config file into `self.cfg`
        Read the yaml into `basedict`
        Store Title, ncpu, constituents, reaction list.
        Create initial_composition list of dict of
            constituent molecules and their counts specified in the YAML.
        Create instances of molecules of initial composition to `self.molecules`
    Set any default values.

Runtime[].do_workflow
    Navigate to proj directory.
    Generate molecules.
        Create instances of all molecules.

    Read the last check point.
        Loads an YAML file and returns dict of sections.
        All do functions have an enableCheckpoint() decorator.

    load_files from the last checkpoint.
        Read the .gro, .top, .grx and .mol2 files into topocoord instance.

    Navigate to systems/init and do_initialization()
        Initialize gromacs topology.
            Take all molecules from the initial composition.
            Charge neutralize the topology.
            Replicate their topology N times, where N is the number of monomers.
            Merge the topology of all molecules.
            Write the 'init.top' file.

        Initialize gromacs coordinates and attributes.
            Take and sanity check the densification yaml configs.
            Calculate boxsize from initial density.
                Take and sanity check the box aspect ratio.
                Calc boxsize using volumen and aspect ratio.
                Set boxsize array in nanometers.

            Let an empty topology dict.
            Take all molecules from the initial composition.
                If the molecule has conformers,
                    Update the topology dict with conformers info until
                    the desired count is reached
                    and copy the files from paremeterized/.
                otherwise, assume a mixture of stereoisomers
                    and copy the files from paremeterized/.

            insert_molecules within boxsize using the topology dict.
                Take the molecule from the inputs/ and
                copy it n times in the box.
            Read the created .gro file.
            Count atoms.
            Update grx attributes using yaml molecules and initial composition.
            Reset id number for cycle and chain.
            Make residue graph.
            Write the grx file.
            Return .gro and .grx file names.

    Navigate to systems/densification and do_densification()
        Perform a series of equilibration specified by the yaml config.
    Navigate to systems/precure and do_precure() and do_cure()
    Navigate to systems/postcure and do_postcure()
    Navigate to systems/final-results and save_data()

```
