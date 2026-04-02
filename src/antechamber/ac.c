int rac(char *filename, int *atomnum, ATOM * atom, int *bondnum,
		BOND * bond, CONTROLINFO *cinfo, MOLINFO *minfo)
{
	int index;
	int tmpint, tmpint1, tmpint2, tmpint3;
	double tmpf1, tmpf2, tmpf3, tmpf4;
	int numatom;
	int numbond;
	int terindex;
	int overflow_flag = 0;
	FILE *fpin;
	char line[MAXCHAR];
	char tmpchar[MAXCHAR];
	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(0);
	}
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	numatom = 0;
	numbond = 0;
	(*cinfo).bpindex = 0;
	terindex = -1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", filename); */
			break;
		}
		if (strncmp("TER", line, 3) == 0) {
			terindex = 1;
			continue;
		}
		if (strncmp("ATOM", line, 4) == 0) {
			if (overflow_flag == 0) {
				if (line[12] != ' ')
					index = 12;
				else
					index = 13;
				atom[numatom].name[0] = line[index];
				atom[numatom].name[1] = line[index + 1];
				atom[numatom].name[2] = line[index + 2];
				atom[numatom].name[3] = line[index + 3];
				atom[numatom].name[4] = '\0';
				index = 17;
				if ((*cinfo).rnindex == 0) {
					atom[numatom].aa[0] = line[index];
					atom[numatom].aa[1] = line[index + 1];
					atom[numatom].aa[2] = line[index + 2];
					atom[numatom].aa[3] = line[index + 3];
					atom[numatom].aa[4] = '\0';
				}
				if (line[21] != ' ')
					atom[numatom].chain[0] = line[21];
				if (terindex == 1) {
					atom[numatom].ter = terindex;
/*the first atom followed after TER line has terindex of 1*/
					terindex = -1;
				}
				sscanf(&line[22], "%d%lf%lf%lf%lf%s",
					   &tmpint, &tmpf1, &tmpf2, &tmpf3, &tmpf4, tmpchar);

				atom[numatom].x = tmpf1;
				atom[numatom].y = tmpf2;
				atom[numatom].z = tmpf3;
				atom[numatom].charge = tmpf4;
				strcpy(atom[numatom].ambername, tmpchar);
			}
			numatom++;

			if (numatom >= (*cinfo).maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
		if (strncmp("BOND", line, 4) == 0) {
			if (overflow_flag == 0) {
				sscanf(&line[4], "%d%d%d%d", &tmpint, &tmpint1,
					   &tmpint2, &tmpint3);
				tmpint1--;
				tmpint2--;
				atom[tmpint1].con[atom[tmpint1].connum++] = tmpint2;
				atom[tmpint2].con[atom[tmpint2].connum++] = tmpint1;
				bond[numbond].bondi = tmpint1;
				bond[numbond].bondj = tmpint2;
				bond[numbond].type = tmpint3;
				if (tmpint3 == 0)
					(*cinfo).bpindex = 1;
			}
			numbond++;
			if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
				printf
					("\nInfo: the bond number exceeds the MAXBOND, reallocate memory, automatically");
				overflow_flag = 1;
			}
		}
		if (strncmp("CHARGE", line, 6) == 0) {
			sscanf(&line[6], "%lf", &(*minfo).dcharge);
			sscanf(&line[18], "%d", &(*minfo).icharge);
		}

	}
	*atomnum = numatom;
	*bondnum = numbond;
/* printf("\n The atomic number is %5d\n", atomnum); */
	fclose(fpin);
	return overflow_flag;

}


/* Antechamber */
void wac(char *filename, int atomnum, ATOM * atom, int bondnum,
		 BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	int i, j, k;
	int num = 0;
	char form[5 * MAXCHAR];

	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		return;
	}
	if (minfo.dcharge < -9990) {
		minfo.dcharge = 0.0;
		for (i = 0; i < atomnum; i++)
			minfo.dcharge += atom[i].charge;
	}
	if (minfo.usercharge < -9990.)
		minfo.icharge = intcharge(atomnum, atom);
	fprintf(fpout, "CHARGE %9.2lf ( %d )\n", minfo.dcharge, minfo.icharge);
	formula(atomnum, atom, form);
	fprintf(fpout, "Formula: %s\n", form);
	for (i = 0; i < atomnum; i++)
		fprintf(fpout,
				"ATOM%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f%10.6lf%10s\n",
				i + 1, atom[i].name, atom[i].aa, atom[i].resno, atom[i].x,
				atom[i].y, atom[i].z, atom[i].charge, atom[i].ambername);

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			for (k = 0; k < bondnum; k++)
				if ((bond[k].bondi == i && bond[k].bondj == j) ||
					(bond[k].bondi == j && bond[k].bondj == i))
					fprintf(fpout, "BOND%5d%5d%5d%5d  %5s%5s\n", 1 + num++,
							i + 1, j + 1, bond[k].type, atom[i].name,
							atom[j].name);

	fclose(fpout);
}
