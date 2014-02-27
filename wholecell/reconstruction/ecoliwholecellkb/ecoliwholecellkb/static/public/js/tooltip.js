/*
 * Tooltip
 * 
 * Author: Jonathan Karr, jkarr@stanford.edu
 * Affiliation: Covert Lab, Department of Bioengineering, Stanford University
 * Last updated: 2012-09-22
 */
 
function showToolTip(evt, title, content){
	tooltip = $('#tooltip');		
	offset = tooltip.parent().offset();	

	tooltip.html('<h1>' + title + '</h1>' + content);
	tooltip.css({left: evt.pageX - tooltip.width() / 2, top: evt.pageY - tooltip.height() - 20})
	if (!tooltip.is(':visible'))
		tooltip.fadeIn(200);
}

function hideToolTip(evt){
	tooltip = $('#tooltip');
	if (tooltip.is(':visible'))
		tooltip.fadeOut(200);
}