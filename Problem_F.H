/* problem-specific stuff goes here */

BL_FORT_PROC_DECL(COM_SAVE,com_save)
    ( const Real* mass_left, const Real* mass_right, 
      const Real* com_xloc_l, const Real* com_xloc_r,
      const Real* com_yloc_l, const Real* com_yloc_r,
      const Real* com_zloc_l, const Real* com_zloc_r
    );

