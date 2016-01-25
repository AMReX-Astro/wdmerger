
#ifdef ROTATION
    //
    // Momemtum
    //
    derive_lst.add("inertial_momentum_x",IndexType::TheCellType(),1,ca_derinertialmomentumx,the_same_box);
    derive_lst.addComponent("inertial_momentum_x",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_momentum_x",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_momentum_y",IndexType::TheCellType(),1,ca_derinertialmomentumy,the_same_box);
    derive_lst.addComponent("inertial_momentum_y",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_momentum_y",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_momentum_z",IndexType::TheCellType(),1,ca_derinertialmomentumz,the_same_box);
    derive_lst.addComponent("inertial_momentum_z",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_momentum_z",desc_lst,State_Type,Xmom,3);

    //
    // Angular momentum
    //
    derive_lst.add("inertial_angular_momentum_x",IndexType::TheCellType(),1,ca_derinertialangmomx,the_same_box);
    derive_lst.addComponent("inertial_angular_momentum_x",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_angular_momentum_x",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_angular_momentum_y",IndexType::TheCellType(),1,ca_derinertialangmomy,the_same_box);
    derive_lst.addComponent("inertial_angular_momentum_y",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_angular_momentum_y",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_angular_momentum_z",IndexType::TheCellType(),1,ca_derinertialangmomz,the_same_box);
    derive_lst.addComponent("inertial_angular_momentum_z",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_angular_momentum_z",desc_lst,State_Type,Xmom,3);
#endif

#ifdef GRAVITY
#ifdef ROTATION
    //
    // Effective potential (gravitational + rotational)
    //
    derive_lst.add("phiEff",IndexType::TheCellType(),1,ca_derphieff,the_same_box);
    derive_lst.addComponent("phiEff",desc_lst,PhiGrav_Type,0,1);
    derive_lst.addComponent("phiEff",desc_lst,PhiRot_Type,0,1);
#endif
#endif
