//replace all document.write with a form-field output or something. document.write is terrible.
if (!String.prototype.format) {
  String.prototype.format = function() {
    var args = arguments;
    return this.replace(/{(\d+)}/g, function(match, number) { 
      return typeof args[number] != 'undefined'
        ? args[number]
        : match
      ;
    });
  };
}

function sequence(header, sequence){
	this.header = header
	this.sequence = sequence
};

function codon(code, letter, mass){
	this.code = code
	this.letter = letter
	this.mass = mass
};

function isolateSubSequences (text){
	var header = /^>.?\|*/i;
	var Seqs= [];
	var lm = 0;
	for (var i=0; i<text.length; i++){
		if (text[i].search(header) != -1){
			console.log("Matched a header at position ", i, ". Last match at ", lm)
			if (lm==0 && i==0){
				continue
			}else{
				Seqs.push(sequence(text[i], text.slice(lm, i-1).join()))
				lm = i
			}			
		};

	};
	if (lm==0){
		Seqs.push(sequence("", text));	
	};
	return Seqs
};


var codons = [
	codon(['AUG'], '@', 0), codon(['UUU', 'UUC'], 'F', 147.06841), codon(['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'L', 113.08406), 
	codon(['AUU', 'AUC', 'AUA'], 'I', 113.08406), codon(['AUG'], 'M', 131.04049), codon(['GUU', 'GUC', 'GUA', 'GUG'], 'V', 99.06841), 
	codon(['UCU', 'UCC', 'UCA', 'UCG'], 'S', 87.03203), codon(['CCU', 'UCC', 'UCA', 'UCG'], 'P', 97.05276), codon(['ACU', 'ACC', 'ACA', 'ACG'], 'T', 101.04768), 
	codon(['GCU', 'GCC', 'GCA', 'GCG'], 'A', 71.0788), codon(['UAU', 'UAC'], 'Y', 163.06333), codon(['UAA', 'UAG', 'UGA'], '*', 0), 
	codon(['CAU', 'CAC'], 'H', 137.05891), codon(['CAA', 'CAG'], 'Q', 128.05858), codon(['AAU', 'AAC'], 'N', 114.04293), 
	codon(['AAA', 'AAG'], 'K', 128.09496), codon(['GAU', 'GAC'], 'D', 115.02694), codon(['GAA', 'GAG'], 'E', 129.04259), 
	codon(['UGU', 'UGC'], 'C', 103.00919), codon(['UGG'], 'W', 186.07931), codon(['AGA', 'AGG', 'CGU', 'CGC', 'CGA', 'CGG'], 'R', 156.10111), 
	codon(['AGU', 'AGC'], 'S', 87.03203), codon(['GGU', 'GGC', 'GGA', 'GGG'], 'G', 57.02146), 
	codon([''], 'B', 114.6068474), codon([''], 'Z', 128.6804964), codon([''], 'X', 110.92618464)];

function translateToProtein (Seq, i){
		//replace T -> U, chop in thirds, case-insensitive match
	protseq=""
	for (var i=0; i<Seq.length; i+3){
		for (var j=0; j<codons.length; j+1){
			if (Seq.slice(i, i+3) == codons[j].code){
				protsec+=codons[j].letter
				console.log("Matched ", Seq.slice(i, i+3), " and ", codons[j].letter)
			};
		};
	};
};


function calculateMass (Seq, i){
	
};

function calculateAbsorbance (Seq, i){
	
};

function parseInput (form) {
	//Parse header/FASTA bullshit
	var header = /^>.?\|*/i;
	//normalize linebreaks and make into an array:
	var rawSeq = form.seqbox.value.replace(/\r\n/g, "\n").split("\n");
	//var Seq = isolateSubSequences(rawSeq)
	if (rawSeq[0].search(header)!=-1){
		
	};
	console.log(rawSeq)
	/*for (var i=0; i<Seq.length; i++){
		Seq[i] = parseSequence(Seq[i], i)
		Seq[i] = calculateMass(Seq[i], i)
		Seq[i] = calculateAbsorbance(Seq[i], i)
	};
	*/
};

function parseSequence (Seq, i){
	var sequenceName = "Unnamed Sequence";
	console.log(Seq, i)
	//Search for non-DNA FASTA characters
	var rxDNA = /[ACGTNUKSYMWRBDHV-]/gi;
	var rxProtein =/[EFILMPQZX*-]/gi;
	var isDNA = Seq.seqeuence.search(rxDNA);
	var isProtein = Seq.search(rxProtein);
	var result;
	var results = ['DNA', 'Protein'];
	switch(isDNA){
		case -1:
			if (isProtein>-1){
				result = 1
				console.log("Is both DNA and protein. Reset to protein.");
			}else{
				result = 0
			};
		default:
			result = 1
	};
		//debug
		document.write("Auto-detecting sequence type for {0}: {1}<br>".format(sequenceName, results[result]));
		//translate if needed
		if (result>0){
			Seq = translateToProtein(Seq);
			};
		//count number of residues
		Res = [];
		for (var i=0; i<Res.length; i++){
			
		};
		//calculate mass
		//calculate coefficient
		//absorbances
			//Wait how do I get absorbances PER sequence
};		
/*
----
var conc = document.form.A.value / ( document.form.e.value * document.form.d.value ) * 1000000;
var mconc = document.form.A.value / ( document.form.e.value * document.form.d.value ) * document.form.M.value;
document.form.concentration.value = conc;
document.form.mconcentration.value = mconc;

----
It has been shown [1c] that it is possible to estimate the molar extinction coefficient of a protein from knowledge of its amino acid composition. From the molar extinction coefficient of tyrosine, tryptophan and cystine (cysteine does not absorb appreciably at wavelengths >260 nm, while cystine does) at a given wavelength, the extinction coefficient of the native protein in water can be computed using the following equation:

E(Prot) = Numb(Tyr)*Ext(Tyr) + Numb(Trp)*Ext(Trp) + Numb(Cystine)*Ext(Cystine)
where (for proteins in water measured at 280 nm): Ext(Tyr) = 1490, Ext(Trp) = 5500, Ext(Cystine) = 125;
The absorbance (optical density) can be calculated using the following formula:
Absorb(Prot) = E(Prot) / Molecular_weight
http://www.ncbi.nlm.nih.gov/pubmed/8563639?dopt=Abstract
----   

	Letter	Amino Acid		Average Mass	% of all residues
	A		alanine			071.0788		8.25
	R		arginine		156.10111		5.53
	N		asparagine		114.04293		4.06
	D		aspartate		115.02694		5.45
	C		cystine			103.00919		1.37
	E		glutamate		129.04259		6.75
	Q		glutamine		128.05858		3.93
	G		glycine			057.02146		7.07
	H		histidine		137.05891		2.27
	I		isoleucine		113.08406		5.95
	L		leucine			113.08406		9.66
	K		lysine			128.09496		5.84
	M		methionine		131.04049		2.42
	F		phenylalanine	147.06841		3.86
	P		proline			097.05276		4.70
	S		serine			087.03203		6.57
	T		threonine		101.04768		5.34
	W		tryptophan		186.07931		1.08
	Y		tyrosine		163.06333		2.92
	V		valine			099.06841		6.87
	U		selenocysteine	150.956363		N/A
	Source: http://web.expasy.org/findmod/findmod_masses.html#AA

	Uncommon codons
	B	aspartate/asparagine 	114,6068474
	Z	glutamate/glutamine		128,6804964
	X	any						110,92618464
						
	These codons are handled as ExPASY's comput_pi (http://web.expasy.org/compute_pi/pi_tool-doc.html) does. Their masses were calculated from ExPASY statistics (http://web.expasy.org/docs/relnotes/relstat.html) on 11th December 2013
	Amino Acid	Letter	Average Mass	% of all residues
	aspartate	D		115.02694		5.45
	asparagine	N		114.04293		4.06
	

	glutamate	E		129.04259		6.75
	glutamine	Q		128.05858		3.93
*/
