# Various scripts to fix MC

## Add lepton masses:
Requires a ROOT installation to run.

This script was written specifically to add lepton masses to the PHANTOM generator for a ZZ decay to electrons/muons. Using the pair of leptons from a Z, a portion of the second lepton is added to the first (and vice versa) such that the resultant 4-momenta have the correct lepton masses and that total momenta of the Z is conserved. _NOTE: Currently, any other generator or process could use this but the code would need considerable tweaking to work properly._

To use, point leptonmassfixDir.sh to the directory containing the LHE files to be altered. This will change ALL LHE files in said directory.

	chmod +x leptonmassfixDir.sh
	./leptonmassfixDir.sh <Directory of LHE files>

