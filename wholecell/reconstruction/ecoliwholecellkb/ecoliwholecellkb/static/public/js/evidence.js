/*
 * Show/hide evidence
 * 
 * Author: Jonathan Karr, jkarr@stanford.edu
 * Affiliation: Covert Lab, Department of Bioengineering, Stanford University
 * Last updated: 2012-09-22
 */
 
function toggleEvidence(self){
	self = $(self);
	
	if (self.parent().find('.evidence').is(':visible')){
		self.parent().find('.evidence').hide(); 
		self.html('(Show evidence)');
	}else{
		self.parent().find('.evidence').show(); 
		self.html('(Hide evidence)');
	}
}