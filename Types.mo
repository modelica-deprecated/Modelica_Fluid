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
end Types;
