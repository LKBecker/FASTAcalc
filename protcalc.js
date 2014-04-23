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

if (!Array.prototype.indexOf) {
    Array.prototype.indexOf = function (obj, fromIndex) {
        if (fromIndex == null) {
            fromIndex = 0;
        } else if (fromIndex < 0) {
            fromIndex = Math.max(0, this.length + fromIndex);
        }
        for (var i = fromIndex, j = this.length; i < j; i++) {
            if (this[i] === obj)
                return i;
        }
        return -1;
    };
}

function commaRound(nnumber, signifFig){
	var snumber = nnumber.toString().split("."), preDecimal = snumber.shift(), postDecimal = snumber.shift();
	var toRound = postDecimal.slice(signifFig-1, postDecimal.length).split("")
	toRound = toRound.map(function(x){return parseInt(x, 10)});
	var finalRounding = toRound.reduceRight(function(last, current){return (current+=(last>4)*1)}, 0);
	if (finalRounding >=10){finalRounding = 9;};
	var rnumber = parseFloat(preDecimal + '.' + postDecimal.slice(0, signifFig-1) + finalRounding);
	return rnumber
};


function sequence(header, sequence){
	this.sequence = String(sequence);
	this.header = String(header);
};

function Codon(codes, letter, name, mass){
	this.codes = codes;
	this.letter = letter;
	this.name = name;
	this.mass = mass;
};

var waterMass = 18.01524;

var codons = [
	new Codon(['GCU', 'GCC', 'GCA', 'GCG'], 				'A', 'Alanine', 71.0788), 
	new Codon(['AGA', 'AGG', 'CGU', 'CGC', 'CGA', 'CGG'], 	'R', 'Arginine', 156.1875), 
	new Codon(['AAU', 'AAC'], 								'N', 'Asparagine', 114.1038), 
	new Codon(['GAU', 'GAC'], 								'D', 'Aspartic acid', 115.0886), 
	new Codon(['UGU', 'UGC'], 								'C', 'Cysteine', 103.1388), 
	new Codon(['GAA', 'GAG'], 								'E', 'Glutamic Acid', 129.1155), 
	new Codon(['CAA', 'CAG'], 								'Q', 'Glutamine', 128.1307), 
	new Codon(['GGU', 'GGC', 'GGA', 'GGG'], 				'G', 'Glycine', 57.0519), 
	new Codon(['CAU', 'CAC'], 								'H', 'Histidine', 137.1411), 
	new Codon(['AUU', 'AUC', 'AUA'], 						'I', 'Isoleucine', 113.1594), 
	new Codon(['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 	'L', 'Leucine', 113.1594), 
	new Codon(['AAA', 'AAG'], 								'K', 'Lysine', 128.1741), 
	new Codon(['AUG'], 										'M', 'Methionine', 131.1926), 
	new Codon(['UUU', 'UUC'], 								'F', 'Phenylalanine', 147.1766), 
	new Codon(['CCU', 'CCC', 'CCA', 'CCG'], 				'P', 'Proline', 97.1167), 
	new Codon(['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 	'S', 'Serine', 87.0782), 
	new Codon(['ACU', 'ACC', 'ACA', 'ACG'], 				'T', 'Threonine', 101.1051), 
	new Codon(['UGG'], 										'W', 'Tryptophan', 186.2132), 
	new Codon(['UAU', 'UAC'], 								'Y', 'Tyrosine', 163.1760), 
	new Codon(['GUU', 'GUC', 'GUA', 'GUG'], 				'V', 'Valine', 99.1326), 
	new Codon(['UAA', 'UAG', 'UGA'], 						'*', 'STOP', 0.0000), 
	new Codon(['XXX'], 										'B', 'Aspartic Acid/Asparagine', 114.6682), 
	new Codon(['SSS'], 										'X', 'Any', 111.3348),
	new Codon(['DDD'], 										'Z', 'Glutamic Acid/Glutamine', 128.7531)];
	
//codons = codons.reduce( function (hash,c){ return c.codes.reduce( function (hash,triplet){ hash[triplet] = c; return hash }, hash)}, {})

function translateProteinSequence (protSeq){
	tokenizedProtSeq = protSeq.split("");
	function letterToCodon(letter){
		matchedLetter = codons.filter( function(c){return (c.letter == letter)} )	
		if (matchedLetter.length >0){
			//console.log(matchedLetter[0]);
			return matchedLetter[0]
		}else{
			return new Codon([],'','',0)
		}
	}
	parsedProtSeq=tokenizedProtSeq.map(letterToCodon)
	return parsedProtSeq
};
	
function translateDNASequence (dnaSeq){
	//translate to RNA
	dnaSeq = dnaSeq.replace(/T/gi, "U")
	dnaSeq = dnaSeq.substring(0, (Math.floor(dnaSeq.length / 3) * 3))
	tokenizedSequence = dnaSeq.split(/(?=(?:...)*$)/)
		
	// function to look up a codon given a string triplet
	function tripletToCodon(triplet) {
		matchedCodons = codons.filter( function(c){return (c.codes.indexOf(triplet) != -1)} )
		if (matchedCodons.length > 0){
			return matchedCodons[0]
		}else{
			// if lookup fails, return a "null" Codon
			return new Codon([],'','',0)
		}
	}
	// turn array of triplets into array of codons
	codonSequence = tokenizedSequence.map(tripletToCodon)
	return codonSequence
};

function calculateExtinctionCoefficient (Seq){
	//for proteins in water measured at 280 nm
	var Extinction = 0
	Tyr = Seq.filter( function(c) {return (c.name=="Tyrosine")});
	Trp = Seq.filter( function(c) {return (c.name=="Tryptophan")});
	Cys = Seq.filter( function(c) {return (c.name=="Cysteine")});
	// console.log(Cys, Trp, Tyr);
	Extinction += (Tyr.length * 1490) + (Trp.length * 5500) + ((Cys.length /2) * 125);
	return Extinction
};

function parseInput (form) {
	//Parse header/FASTA bullshit
	var header = /^>.?\|*/i;
	var sheader = "(Unnamed Sequence)"
	//normalize linebreaks and make into an array:
	var rawSeq = form.seqbox.value.replace(/\r\n/g, "\n").split("\n");
	//var Seq = isolateSubSequences(rawSeq)
	if (rawSeq[0].search(header)!=-1){
		sheader = rawSeq.slice(0, 1).join("");
		rawSeq = rawSeq.slice(1, rawSeq.length+1);
	};
	Seq = new sequence(sheader, rawSeq.join(""));
	Seq.rawSequence = Seq.sequence;
	var rawSeqType;
	if (form.FASTAtype[0].checked){
		Seq.sequence = translateProteinSequence(Seq.sequence);
	}else{
		Seq.sequence = translateDNASequence(Seq.sequence);
	}
	Seq.strSequence = Seq.sequence.map( function(c){ return c.letter} ).join("");
	Seq.mass = Seq.sequence.reduce( function(totalMass, c){ return totalMass + c.mass}, 0);
	Seq.mass += waterMass
	Seq.extinction = calculateExtinctionCoefficient(Seq.sequence);
	Seq.absorbance = parseFloat(form.Absorbance.value);
	Seq.pathlength = parseFloat(form.l.value);
	// console.log(Seq.header, Seq.mass, Seq.extinction, Seq.absorbance, Seq.pathlength);
	
	Seq.molConc = (Seq.absorbance/(Seq.pathlength * Seq.extinction))
	
	var sigFig = parseInt(form.signFig.value);
	// console.log(sigFig, Seq.extinction, commaRound(Seq.extinction, sigFig), Seq.extinction.toFixed(sigFig));
	
	
	if (form.FASTAtype[1].checked){
		document.getElementById("protSeq").innerHTML=Seq.strSequence;
		document.getElementById("translation").style.display="table"};
	
	document.getElementById("protName").innerHTML=Seq.header;
	document.getElementById("molMass").innerHTML=Seq.mass.toFixed(sigFig);
	document.getElementById("extCoeff").innerHTML=Seq.extinction.toFixed(sigFig);
	document.getElementById("umolConc").innerHTML=(Seq.molConc * 1000000).toFixed(sigFig);
	document.getElementById("mgMlConc").innerHTML=(Seq.molConc * Seq.mass).toFixed(sigFig);
	document.getElementById("results").style.display="table";
};

function blank (form){
	form.seqbox.value = "";
	form.Absorbance.value = "";
	form.l.value = 1;
	form.signFig.value = 1;
	form.FASTAtype[0].checked = true;
	document.getElementById("results").style.display="none";
	document.getElementById("translation").style.display="none";
};

/*
----
var conc = document.form.A.value / ( document.form.e.value * document.form.d.value ) * 1000000;
var mconc = document.form.A.value / ( document.form.e.value * document.form.d.value ) * document.form.M.value;
document.form.concentration.value = conc;
document.form.mconcentration.value = mconc;

----
It has been shown that it is possible to estimate the molar extinction coefficient of a protein from knowledge of its amino acid composition. From the molar extinction coefficient of tyrosine, tryptophan and cystine (cysteine does not absorb appreciably at wavelengths >260 nm, while cystine does) at a given wavelength, the extinction coefficient of the native protein in water can be computed using the following equation:

E(Prot) = Numb(Tyr)*Ext(Tyr) + Numb(Trp)*Ext(Trp) + Numb(Cystine)*Ext(Cystine)
where (for proteins in water measured at 280 nm): Ext(Tyr) = 1490, Ext(Trp) = 5500, Ext(Cystine) = 125;

The absorbance (optical density) can be calculated using the following formula:
Absorb(Prot) = E(Prot) / Molecular_weight
http://www.ncbi.nlm.nih.gov/pubmed/8563639?dopt=Abstract
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2143013/
----   

	Letter	Amino Acid		Average mass	% of all residues
	A		alanine			071.0788		8.25%
	R		arginine		156.1875		5.53%
	N		asparagine		114.1038		4.06%
	D		aspartate		115.0886		5.45%
	C		cystine			103.1388		1.37%
	E		glutamate		129.1155		6.75%
	Q		glutamine		128.1307		3.93%
	G		glycine			057.0519		7.07%
	H		histidine		137.1411		2.27%
	I		isoleucine		113.1594		5.95%
	L		leucine			113.1594		9.66%
	K		lysine			128.1741		5.84%
	M		methionine		131.1926		2.41%
	F		phenylalanine	147.1766		3.86%
	P		proline			097.1167		4.70%
	S		serine			087.0782		6.57%
	T		threonine		101.1051		5.34%
	W		tryptophan		186.2132		1.08%
	Y		tyrosine		163.176			2.92%
	V		valine			099.1326		6.86%
	U		selenocysteine	150.0388		N/A	
	Source: http://web.expasy.org/findmod/findmod_masses.html#AA

	Uncommon amino acids
	B	aspartate/asparagine 	114.6682
	Z	glutamate/glutamine		128.7531
	X	any						110.9698
						
	These are handled as ExPASY's comput_pi (http://web.expasy.org/compute_pi/pi_tool-doc.html) does. Their masses were calculated from ExPASY statistics (http://web.expasy.org/docs/relnotes/relstat.html) on 12th March 2014.
	Amino Acid	Letter	Average Mass	% of all residues
	aspartate	D		115.0886		5.45
	asparagine	N		114.1038		4.06
	
	glutamate	E		129.1155		6.75
	glutamine	Q		128.1307		3.93
*/
