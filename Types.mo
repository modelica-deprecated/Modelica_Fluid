package Types "Common types for fluid models" 
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
Package <b>Types</b> contains common type definitions of the Modelica_Fluid
library.
</p>
</html>"));
  
  type DensityUnits "Menu choices for density units" 
      extends Modelica.Icons.TypeString;
      annotation (Evaluate=true,
                 choices(choice="kg/m3" "\"kg/m3\"",
                         choice="g/cm3" "\"g/cm3\" (= 1e-3 kg/m3)"));
  end DensityUnits;
  
  type MassFlowRateUnits "Menu choices for mass flow rate units" 
      extends Modelica.Icons.TypeString;
      annotation (Evaluate=true,
                 choices(choice="kg/s" "\"kg/s\"",
                         choice="t/h" "\"t/h\""));
  end MassFlowRateUnits;
  
  type PressureUnits "Menu choices for pressure units" 
      extends Modelica.Icons.TypeString;
      annotation (Evaluate=true,
                 choices(choice="Pa" "\"Pa\" (Pascal)",
                         choice="kPa" "\"kPa\" (= 1e3 Pa)",
                         choice="bar" "\"bar\" (= 1e5 Pa)",
                         choice="MPa" "\"Mpa\" (= 1e6 Pa)"));
    
  end PressureUnits;
  
  type SpecificEnthalpyUnits "Menu choices for specific enthalpy units" 
      extends Modelica.Icons.TypeString;
      annotation (Evaluate=true,
                 choices(choice="J/kg" "\"J/kg\"",
                         choice="kJ/kg" "\"kJ/kg\" (= 1e3 J/kg)"));
  end SpecificEnthalpyUnits;
  
  type TemperatureUnits "Menu choices for temperature units" 
      extends Modelica.Icons.TypeString;
      annotation (Evaluate=true,
                 choices(choice="K" "\"K\" (Kelvin)",
                         choice="degC" "\"degC\" (degree Celsius)",
                         choice="degF" "\"degF\" (degree Fahrenheit)",
                         choice="degR" "\"degR\" (degree Rankine)"));
  end TemperatureUnits;
  
  type VolumeFlowRateUnits "Menu choices for volume flow rate units" 
      extends Modelica.Icons.TypeString;
      annotation (Evaluate=true,
                 choices(choice="m3/s" "\"m3/s\"",
                         choice="l/s" "\"l/s\""));
  end VolumeFlowRateUnits;
  
  package FrictionTypes 
    "Type, constants and menu choices to define the pressure loss equations due to friction, as temporary solution until enumerations are available" 
    
    extends Modelica.Icons.Library;
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
    
    extends Modelica.Icons.Library;
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
    
    extends Modelica.Icons.Library;
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
end Types;
