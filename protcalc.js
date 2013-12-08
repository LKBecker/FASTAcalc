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

function isolateSubsequence (text){
	Seqs = [];
	for (var i=0; i<text.length; i++){
		if (text[i].search(header) != -1){
			console.log("Matched a header at position ", i, ". Last match at ", lm)
			Seqs.push(text.slice(lm, i-1).join())
			lastmatch = i
		};
	};
	return Seqs
};

function translateToProtein (text){
	var codons=[];
	for (var i=0; i<text.length; i+3){
	
	};
};
	
function parseInput (form) {
	//Parse header/FASTA bullshit
	//normalize linebreaks and make into an array:
	var header = /^>.?\|*/i;
	var rawSeqs = form.seqbox.value.replace(/\r\n/g, "\n").split("\n");
	var Seqs = [];
	var lm = 0;

	var sequenceName = "Unnamed Sequence";
	console.log("Seqs: ", Seqs);
	document.write("{0} sequence(s) found.<br>".format(Seqs.length));
	for (var i=0; i<Seqs.length; i++){
		Seq = Seqs[i]	
		//Search for non-DNA FASTA characters
		var rxDNA = /[ACGTNUKSYMWRBDHV-]/gi;
		var rxProtein =/[EFILMPQZX*-]/gi;
		var isDNA = Seq.search(rxDNA);
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
		//count number of residues
		//calculate mass
		//calculate coefficient
		//absorbances
			//Wait how do I get absorbances PER sequence
	};
		
};

	

