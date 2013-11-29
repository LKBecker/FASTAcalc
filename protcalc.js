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

function translateToProtein (text) {

	};
	
function parseInput (form) {
	//Parse header/FASTA bullshit
	//normalize linebreaks and make into an array:
	var header = /^>.?\|*/i;
	var rawSeqs = form.seqbox.value.replace(/\r\n/g, "\n").split("\n");
	var Seqs = [];
	var lm = 0;
	for (var i = 0; i<rawSeqs.length; i++){
		if (rawSeqs[i].search(header) != -1){
			console.log("Matched a header at position ", i, ". Last match at ", lm)
			Seqs.push(rawSeqs.slice(lm, i-1).join())
			lastmatch = i
		}
	}
	var sequenceName = "Unnamed Sequence";
	console.log("Seqs: ", Seqs);
	document.write("{0} sequences found.<br>".format(Seqs.length));
	for (var i = 0; i<Seqs.length; i++)
		{
			
			
		document.write("Auto-detecting sequence type for {0}:<br>".format(sequenceName));
		//Search for non-DNA FASTA characters
		var rxDNA = /[ACGTNUKSYMWRBDHV-]/gi;
		var rxProtein =/[EFILMPQZX*-]/gi;
		var isDNA = Seq.search(rxDNA);
		var isProtein = Seq.search(rxProtein);
		if (isDNA>-1 && isProtein>-1){
			console.log("Is both DNA and protein. Reset to protein.");
			isDNA = -1;
			isProtein = 1;
		};
		
		//debug
		document.write("<t>Sequence is DNA: {0}<br>".format(isDNA>-1));
		document.write("<t>Sequence is Protein: {0}<br>".format(isProtein>-1));
		};
	};
	

