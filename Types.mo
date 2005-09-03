package Types "Common types for fluid models" 
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
Package <b>Types</b> contains common type definitions of the Modelica_Fluid
library.
</p>
</html>"));
  
  package FrictionTypes 
    "Type, constants and menu choices to define the pressure loss equations due to friction, as temporary solution until enumerations are available" 
    
    extends Modelica.Icons.Enumeration;
    constant Integer ConstantLaminar=1;
    constant Integer ConstantTurbulent=2;
    constant Integer DetailedFriction=3;
    type Temp 
      "Temporary type of FrictionTypes with choices for menus (until enumerations are available)" 
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar 
            "ConstantLaminar \"dp = k*m_flow\"",
          choice=Modelica_Fluid.Types.FrictionTypes.ConstantTurbulent 
            "ConstantTurbulent \"dp = k*m_flow^2\"",
          choice=Modelica_Fluid.Types.FrictionTypes.DetailedFriction 
            "DetailedFriction \"dp = f(Re,delta,rho,L,D,nu)\""));
    end Temp;
  end FrictionTypes;
  
  package CrossSectionTypes 
    "Type, constants and menu choices to define the geometric cross section of pipes, as temporary solution until enumerations are available" 
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Enumeration;
    constant Integer Circular=1;
    constant Integer Rectangular=2;
    constant Integer General=3;
    type Temp 
      "Temporary type of CrossSectionTypes with choices for menus (until enumerations are available)" 
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.CrossSectionTypes.Circular 
            "Circular cross section",
          choice=Modelica_Fluid.Types.CrossSectionTypes.Rectangular 
            "Rectangular cross section",
          choice=Modelica_Fluid.Types.CrossSectionTypes.General 
            "General cross section"));
    end Temp;
  end CrossSectionTypes;
  
  package InitTypes 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit=1;
    constant Integer InitialStates=2;
    constant Integer SteadyState=3;
    constant Integer SteadyMass=4;
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.InitTypes.NoInit 
            "NoInit (no initialization)",
          choice=Modelica_Fluid.Types.InitTypes.InitialStates 
            "InitialStates (initialize medium states)",
          choice=Modelica_Fluid.Types.InitTypes.SteadyState 
            "SteadyState (initialize in steady state)",
          choice=Modelica_Fluid.Types.InitTypes.SteadyMass 
            "SteadyMass (initialize density or pressure in steady state)"));
    end Temp;
  end InitTypes;
  
  model CvTypes 
    "Type, constants and menu choices to define the choice of valve flow coefficient" 
    extends Modelica.Icons.Enumeration;
    annotation (preferedView="text");
    constant Integer Av = 0 "Av (metric) flow coefficient";
    constant Integer Kv = 1 "Kv (metric) flow coefficient";
    constant Integer Cv = 2 "Cv (US) flow coefficient";
    constant Integer OpPoint = 3 "Av defined by initial equation";
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=0, max=3);
      annotation (Evaluate=true, choices(
        choice=Modelica_Fluid.Types.CvTypes.Av "Av (metric) flow coefficient",
        choice=Modelica_Fluid.Types.CvTypes.Kv "Kv (metric) flow coefficient",
        choice=Modelica_Fluid.Types.CvTypes.Cv "Cv (US) flow coefficient",
        choice=Modelica_Fluid.Types.CvTypes.OpPoint 
            "Av defined by nominal operating point"));
    end Temp;
    
  end CvTypes;
end Types;
