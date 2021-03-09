 image=getTitle;
print(image);
setTool("multipoint");
while(1) {
	makePoint(0, 0);
	run("Properties... ", "  stroke=red point=Circle size=[Extra Large] show");
	waitForUser("Click OK to clear");
	selectWindow("Counts_"+image);
	info = split(getInfo(),'\n');
	print(substring(info[1],5,lengthOf(info[1])));
}