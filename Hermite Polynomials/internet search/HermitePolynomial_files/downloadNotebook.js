function downloadNotebook(url) {
	if(readCookie('gotmathematica') != null) {
		document.location = url;
	} else {
		window.open('http://demonstrations.wolfram.com/player.html?url=' + url,'PlayerPopup','toolbar=no,location=no,directories=no,status=no,menubar=no,scrollbars=no,resizable=no,width=456,height=284,left=50,top=50');
		return false;	
	}
}

function readCookie(name) {
	var nameEQ = name + "=";
	var ca = document.cookie.split(';');
	for(var i=0;i < ca.length;i++) {
		var c = ca[i];
		while (c.charAt(0)==' ') c = c.substring(1,c.length);
		if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
	}
	return null;
}