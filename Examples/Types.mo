package Types 
  import SI = Modelica.SIunits;
  annotation (preferedView="text");
  
  type Length_mm = Real (
      quantity="Length", 
      unit="mm", 
      min=0);
  
  
    // constant SI.MassFlowRate ResidualFlow=1E-10 "Used to take care of zero mass flow";
  
  package FrictionTypes 
    "Type, constants and menu choices to define the pressure loss equations due to friction, as temporary solution until enumerations are available"
     
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer ConstantLaminar=1;
    constant Integer ConstantTurbulent=2;
    constant Integer DetailedFriction=3;
    type Temp 
      "Temporary type of FrictionTypes with choices for menus (until enumerations are available)"
       
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar 
            "ConstantLaminar \"dp = k*m_dot\"", 
          choice=Modelica_Fluid.Types.FrictionTypes.ConstantTurbulent 
            "ConstantTurbulent \"dp = k*m_dot^2\"", 
          choice=Modelica_Fluid.Types.FrictionTypes.DetailedFriction 
            "DetailedFriction \"dp = f(Re,delta,rho,L,D,nu)\""));
    end Temp;
  end FrictionTypes;
  
  package CrossSectionTypes 
    "Type, constants and menu choices to define the geometric cross of pipes, as temporary solution until enumerations are available"
     
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer Circular=1;
    constant Integer Rectangular=2;
    constant Integer General=3;
    type Temp 
      "Temporary type of CrossSectionTypes with choices for menus (until enumerations are available)"
       
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.CrossSectionTypes.Circular "Circular", 
          choice=Modelica_Fluid.Types.CrossSectionTypes.Rectangular 
            "Rectangular", 
          choice=Modelica_Fluid.Types.CrossSectionTypes.General "General"));
    end Temp;
  end CrossSectionTypes;
  
  package InitTypes 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available"
     
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer SteadyState=1;
    constant Integer SteadyPressure=2;
    constant Integer InitialStates=3;
    constant Integer NoDefaultInit=4;
    type Temp 
      "Temporary type of InitializationTypes with choices for menus (until enumerations are available)"
       
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.InitializationTypes.SteadyState 
            "SteadyState (initialize in steady state)", 
          choice=Modelica_Fluid.Types.InitializationTypes.SteadyPressure 
            "SteadyPressure (initialize pressure in steady state)", 
          choice=Modelica_Fluid.Types.InitializationTypes.InitialStates 
            "InitialStates (initialize medium states)", 
          choice=Modelica_Fluid.Types.InitializationTypes.NoDefaultInit 
            "NoDefaultInit (no default initialization)"));
    end Temp;
  end InitTypes;
  
  package FixedDensityTypes 
    "Type, constants and menu choices to define fixed density definition, as temporary solution until enumerations are available"
     
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer NoFixedDensity=1;
    constant Integer FixedDensity=2;
    constant Integer FixedDensity_via_pT=3;
    type Temp 
      "Temporary type of FixedDensityTypes with choices for menus (until enumerations are available)"
       
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.FixedDensityTypes.NoFixedDensity 
            "NoFixedDensity (no fixed density for compressible media)", 
          choice=Modelica_Fluid.Types.FixedDensityTypes.FixedDensity 
            "FixedDensity (fixed density for compressible media)", 
          choice=Modelica_Fluid.Types.FixedDensityTypes.FixedDensity_via_pT 
            "FixedDensity_via_pT (fixed density defined by p, T)"));
    end Temp;
  end FixedDensityTypes;
  
end Types;
