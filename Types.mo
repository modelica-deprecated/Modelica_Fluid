package Types "Common types for fluid models" 
  type TemperatureUnits "Menu choices for temperature units" 
      extends String;
      annotation (Evaluate=true,
                 choices(choice="K" "\"K\" (Kelvin)",
                         choice="degC" "\"degC\" (degree Celsius)",
                         choice="degF" "\"degF\" (degree Fahrenheit)",
                         choice="degR" "\"degR\" (degree Rankine)"), 
      Icon(Rectangle(extent=[-100,100; 100,-100], style(
            color=0, 
            rgbcolor={0,0,0}, 
            fillColor=8, 
            rgbfillColor={181,181,181})), Text(
          extent=[-100,100; 100,-100], 
          string="T", 
          style(
            color=0, 
            rgbcolor={0,0,0}, 
            fillColor=8, 
            rgbfillColor={181,181,181}))));
  end TemperatureUnits;
  
  package PressureUnits 
    "Type, constants and menu choices for pressure units, as temporary solution until enumerations will be available" 
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer Pascal=1;
    constant Integer Bar=2;
    constant Integer KiloPascal=3;
    constant Integer MegaPascal=4;
    type Temp 
      "Temporary type of PressureUnits with choices for menus (until enumerations will be available)" 
      
      extends Integer;
      annotation (Evaluate=true,
                 choices(choice=Modelica_Fluid.Types.PressureUnits.Pascal "Pa",
                         choice=Modelica_Fluid.Types.PressureUnits.Bar "bar",
                         choice=Modelica_Fluid.Types.PressureUnits.KiloPascal "kPa",
                         choice=Modelica_Fluid.Types.PressureUnits.MegaPascal "MPa"));
    end Temp;
  end PressureUnits;
  import SI = Modelica.SIunits;
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
This package will contain common type definitions of the fluid library.
It is currently empty as all types are being evaluated together 
with examples first (see sub-package Examples). 
</p>
</html>"));
  
  package MassFlowRateUnits 
    "Type, constants and menu choices for mass flow rate units, as temporary solution until enumerations will be available" 
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer KilogramPerSecond=1;
    constant Integer TonPerHour=2;
    type Temp 
      "Temporary type of MassFlowRateUnits with choices for menus (until enumerations will be available)" 
      
      extends Integer;
      annotation (Evaluate=true,
                 choices(choice=Modelica_Fluid.Types.MassFlowRateUnits.KilogramPerSecond "kg/s",
                         choice=Modelica_Fluid.Types.MassFlowRateUnits.TonPerHour "t/h"));
    end Temp;
  end MassFlowRateUnits;
end Types;
