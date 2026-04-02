/* PDB */
int rpdb(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		 MOLINFO minfo, int pqr)
{
	int numatom;
	int tmpint2;
	int terindex;
	int overflow_flag = 0;
	char tmpchar[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1, tmpfloat2;
	double x, y, z;
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	numatom = 0;
	initial(cinfo.maxatom, atom, minfo.resname);
	terindex = -1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
			/*  printf("\nFinished reading %s file.", filename); */
			break;
		}
		if (strncmp("TER", line, 3) == 0) {
			terindex = 1;
			continue;
		}
		if (strncmp("ATOM", line, 4) == 0
			|| strncmp("HETATM", line, 6) == 0) {
			if (overflow_flag == 0) {
				if (pqr)
					sscanf(&line[22], "%d%lf%lf%lf%lf%lf%s", &tmpint2, &x,
						   &y, &z, &tmpfloat1, &tmpfloat2, tmpchar);
				else
					sscanf(&line[22], "%d%lf%lf%lf", &tmpint2, &x, &y, &z);

/*          --- columns 13-16 have the Brookhaven-formatted name:    */

				atom[numatom].name[0] = line[12];
				atom[numatom].name[1] = line[13];
				atom[numatom].name[2] = line[14];
				atom[numatom].name[3] = line[15];

/*          --- now unwrap this to a more understandable convention:    */

				if (atom[numatom].name[0] == ' ' ||
					atom[numatom].name[0] == '1' ||
					atom[numatom].name[0] == '2' ||
					atom[numatom].name[0] == '3') {

					atom[numatom].name[0] = line[13];
					atom[numatom].name[1] = line[14];
					atom[numatom].name[2] = line[15];
					atom[numatom].name[3] = line[12];

				}
				atom[numatom].name[4] = '\0';
				if (cinfo.rnindex == 0) {
					atom[numatom].aa[0] = line[17];
					atom[numatom].aa[1] = line[18];
					atom[numatom].aa[2] = line[19];
					atom[numatom].aa[3] = '\0';
				}
				if (line[21] != ' ')
					atom[numatom].chain[0] = line[21];
				atom[numatom].ter = terindex;
				atom[numatom].resno = tmpint2;
				atom[numatom].x = x;
				atom[numatom].y = y;
				atom[numatom].z = z;
				if (terindex == 1) {
					atom[numatom].ter = terindex;
					terindex = -1;
				}
				if (pqr) {
					atom[numatom].charge = tmpfloat1;
					atom[numatom].radius = tmpfloat2;
					strcpy(atom[numatom].ambername, tmpchar);
				}

				if (strcmp(atom[numatom].name, "dumm") == 0)
					continue;
				if (strcmp(atom[numatom].name, "Du") == 0)
					continue;
				if (strcmp(atom[numatom].name, "DUMM") == 0)
					continue;
			}
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nThe atom number exceeds the MAXATOM, reallocate memory");
				overflow_flag = 1;
			}
		}
	}
	*atomnum = numatom;
	fclose(fpin);
	return overflow_flag;
}



void wpdb(char *filename, int atomnum, ATOM atom[])
{
	FILE *fpout;
	int i;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atomnum; i++)
		fprintf(fpout, "ATOM%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f\n",
				i + 1, atom[i].name, atom[i].aa, atom[i].resno, atom[i].x,
				atom[i].y, atom[i].z);
	fclose(fpout);
}

void wmpdb(char *filename, int atomnum, ATOM atom[])
{
	FILE *fpout;
	int i;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atomnum; i++)
		fprintf(fpout,
				"ATOM%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f%10.6lf%8.2lf%8s\n",
				i + 1, atom[i].name, atom[i].aa, atom[i].resno, atom[i].x,
				atom[i].y, atom[i].z, atom[i].charge, atom[i].radius,
				atom[i].ambername);
	fclose(fpout);
}
