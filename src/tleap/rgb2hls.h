#ifndef __RGBtoHLS__H
#define __RGBtoHLS__H

#if NeedFunctionPrototypes
#ifndef _
#define _(args) args
#endif
#else
#ifndef _
#define _(args) ()
#endif
#endif

extern void RGBtoHLS _( (double,
			 double,
			 double,
			 double*,
			 double*,
			 double*)
		       );
		       
extern void HLStoRGB _( (double*,
			 double*,
			 double*,
			 double,
			 double,
			 double)
		       );

#ifdef NeedShadowColor

extern Pixel  TopShadowColor(
#if NeedFunctionPrototypes
 Widget  /* self */,
 Pixel   /* base */
#endif
);

extern Pixel  BottomShadowColor(
#if NeedFunctionPrototypes
 Widget  /* self */,
 Pixel   /* base */
#endif
);

#endif /* NeedShadowColor */
#endif /* __RGBtoHLS__H */
