void memory(int flag, int maxatom, int maxbond, int maxring)
{
	if (flag == 0) {
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "meory allocation error for *atom\n");
			exit(0);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "meory allocation error for *arom\n");
			exit(0);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "meory allocation error for *bond\n");
			exit(0);
		}
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stderr, "meory allocation error for *ring\n");
			exit(0);
		}
	}
	if (flag == 3) {
		free(atom);
		free(arom);
		free(bond);
		free(ring);
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "meory allocation error for *atom\n");
			exit(0);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "meory allocation error for *arom\n");
			exit(0);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "meory allocation error for *bond\n");
			exit(0);
		}
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stderr, "meory allocation error for *ring\n");
			exit(0);
		}
	}
	if (flag == 1) {
		free(atom);
		free(arom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stderr, "meory allocation error for *atom\n");
			exit(0);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stderr, "meory allocation error for *arom\n");
			exit(0);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stderr, "meory allocation error for *bond\n");
			exit(0);
		}
	}
	if (flag == 2) {
		free(ring);
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stderr, "meory allocation error for *ring\n");
			exit(0);
		}
	}
}
