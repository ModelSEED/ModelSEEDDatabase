import java.io.IOException;
import java.util.HashMap;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;

import chemaxon.standardizer.Standardizer;
import chemaxon.standardizer.configuration.StandardizerConfiguration;
import chemaxon.standardizer.actions.RemoveExplicitHydrogensAction;
import chemaxon.standardizer.actions.AromatizeAction;
import chemaxon.standardizer.actions.ConvertPiMetalBondsAction;
import chemaxon.standardizer.actions.DisconnectMetalAtomsAction;
import chemaxon.standardizer.actions.TautomerizeAction;
import chemaxon.standardizer.actions.SetAbsoluteStereoAction;

public final class StandardizeMol {

    public static void main(String[] args){

	// Initialize Standardizer Configuration with set of actions
	// This list was partly based on PubChem:
	// https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0293-8

	// All action parameters are empty because there is no documentation on the parameters
	HashMap<String, String> parameters = new HashMap<String, String>();
	StandardizerConfiguration standardizer_configuration = new StandardizerConfiguration();
	standardizer_configuration.addAction(new RemoveExplicitHydrogensAction(parameters));
	standardizer_configuration.addAction(new ConvertPiMetalBondsAction(parameters));
	standardizer_configuration.addAction(new DisconnectMetalAtomsAction(parameters));
	standardizer_configuration.addAction(new AromatizeAction(parameters));
	standardizer_configuration.addAction(new TautomerizeAction(parameters));
	standardizer_configuration.addAction(new SetAbsoluteStereoAction(parameters));

	// Initialize Standardizer
	Standardizer standardizer = new Standardizer(standardizer_configuration);
		
	//Import Mol File
	Molecule OriginalMol = new Molecule();
	try {
	    
	    MolImporter mi = new MolImporter(args[0]);
	    
	    try {
		OriginalMol = mi.read();
	    } finally {
		mi.close();
	    }
	    
	} catch (MolFormatException e) {
	    throw new IllegalArgumentException("Invalid molecule format", e);
	} catch (IOException e) {
	    throw new IllegalArgumentException("Error reading input file", e);
	}
	
	//Molecule for fusing all fragments
	Molecule FusedMol = new Molecule();

	//Fragment Molecule and iterate
	Molecule frags[] = OriginalMol.convertToFrags();
	for(int i = 0; i < frags.length; i++){
	    
	    standardizer.standardize(frags[i]);
	    
	    frags[i].dearomatize();
	    
	    //Fuse fragment
	    FusedMol.fuse(frags[i],false);

	}

	//Print to string
	try {
	    String FusedString = MolExporter.exportToFormat(FusedMol, "mol");
	    System.out.print(FusedString);

	} catch (IOException IOE) {
	    System.out.println("Error: "+IOE.getMessage());
	}
    }
}
