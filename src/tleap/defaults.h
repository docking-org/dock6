#define MAXRESID 50000

typedef struct {
	double	dDSearchDistance;
	double	dDESPGridSpace;
	double	dDESPBoxSize;
	float	fDESPDielectric;
	int	iDESPConstant;
	int	pdbwritecharges;
	double	dGridSpace;
	double	dShellExtent;
	int	iDielectricFlag;
	int	iGBparm;
	int 	iOldPrmtopFormat;
	int		iGibbs;
	int 	iCharmm; 
	int	iResidueImpropers;
	int	iDeleteExtraPointAngles;
	int	iPdbLoadSequential;
	int	iHaveResIds;
	int iUseResIds;
	int iFlexibleWater;
	char	sResidueId[MAXRESID][7];
} defaultstruct ;

extern defaultstruct GDefaults; 

