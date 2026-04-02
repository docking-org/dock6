/* GCRT */
int rgcrt(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		  MOLINFO minfo)
{
	FILE *fpin;
	int i, index;
	int nameindex;
	int numatom;
	int overflow_flag = 0;
	char tmpchar[10];
	char line[MAXCHAR];


	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	index = 0;
	numatom = 0;
	nameindex = -1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*     printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		if (spaceline(line) == 1) {
			index++;
			continue;
		}
		if (index >= 2)
			index++;
		if (index <= 3)
			continue;
		if (index > 4 && spaceline(line) == 1)
			break;
		sscanf(line, "%s", tmpchar);
		if (nameindex == -1) {
			if (tmpchar[0] >= '0' && tmpchar[0] <= '9')
				nameindex = 0;
			else
				nameindex = 1;
		}
		if (overflow_flag == 0) {
			if (nameindex == 0)
				sscanf(line, "%d%lf%lf%lf", &atom[numatom].atomicnum,
					   &atom[numatom].x, &atom[numatom].y,
					   &atom[numatom].z);
			else
				sscanf(line, "%s%lf%lf%lf", atom[numatom].name,
					   &atom[numatom].x, &atom[numatom].y,
					   &atom[numatom].z);
		}
		numatom++;
		if (numatom >= cinfo.maxatom && overflow_flag == 0) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
	}
	*atomnum = numatom;
	if (nameindex == 0) {
		element(*atomnum, atom);
		for (i = 0; i < *atomnum; i++)
			strcpy(atom[i].name, atom[i].element);
	}
	if (nameindex == 1)
		atomicnum(*atomnum, atom);
/* printf("\n atom number is  %5d", *atomnum); */
	fclose(fpin);
	return overflow_flag;
}

void wgcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i;
	/* int index; */

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	fprintf(fpout, "%s\n", "--Link1--");
	fprintf(fpout, "%s%s\n", "%chk=", minfo.chkfile);
	fprintf(fpout, "%s\n\n", minfo.gkeyword);
	fprintf(fpout, "%s\n\n", "remark line goes here");
	fprintf(fpout, "%d%4d\n", minfo.icharge, minfo.multiplicity);
	element(atomnum, atom);

	for (i = 0; i < atomnum; i++)
		fprintf(fpout, "%5s%12.4lf    %12.4lf    %12.4lf     \n",
				atom[i].element, atom[i].x, atom[i].y, atom[i].z);

	fprintf(fpout, "\n");
	fclose(fpout);
}
