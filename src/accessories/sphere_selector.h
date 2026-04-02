using namespace std;

typedef struct atom_def{
	int		idx;
	float	x;
	float	y;
	float	z;
	float	radius;
	int		idx2;
	bool	selected;
}	SPHERE;

typedef	struct coords_def{
	float	coord[3];
}	XYZ;

class Spheres{
	
public:

		~Spheres();

		vector< SPHERE > tot_spheres;
		int			num_spheres, atoms, num_selected;
		XYZ			*atom_coord;

		void		read_spheres(FILE *);
		void		read_atoms(FILE *);
		float		distance(int, int);
		void		print_output(char *);

};
