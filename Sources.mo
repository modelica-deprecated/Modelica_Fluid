package Sources "Generic fluid sources" 
  extends Modelica.Icons.Library;
  model FixedAmbient "Ambient temperature and pressure source" 
    import Modelica.SIunits.Conversions.*;
    
    extends Interfaces.PartialSource(medium(known_pd=if use_p_ambient or Medium.
             incompressible then Medium.Choices.pd.p_known else Medium.Choices.
            pd.d_known, known_Th=if use_T_ambient then Medium.Choices.Th.
            T_known else Medium.Choices.Th.h_known,
            init_p = true,
            p_start = p_ambient,
            T_start = T_ambient,
            X_start = X_ambient,
            h_start = h_ambient));
    
    parameter Boolean use_p_ambient=true 
      "|Ambient pressure or ambient density| = true, if p_ambient is used, otherwise d_ambient (true is required for incompressible medium)"
      annotation (Evaluate=true);
    parameter Modelica_Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=101325 " Ambient pressure, if use_p_ambient = true" annotation 
      (Dialog(group="Ambient pressure or ambient density", enable=use_p_ambient));
    parameter Modelica_Media.Interfaces.PartialMedium.Density d_ambient=1 
      " Ambient density, if use_p_ambient = false" annotation (Dialog(group=
            "Ambient pressure or ambient density", enable=not use_p_ambient));
    parameter Boolean use_T_ambient=true 
      "|Ambient temperature or ambient specific enthalpy| = true, if T_ambient is used, otherwise h_ambient"
      annotation (Evaluate=true);
    parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=
        from_degC(20) " Ambient temperature, if use_T_ambient = true" 
      annotation (Dialog(group=
            "Ambient temperature or ambient specific enthalpy", enable=
            use_T_ambient));
    parameter Modelica_Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=1.e4 " Ambient specific enthalpy, if use_T_ambient = false" 
      annotation (Dialog(group=
            "Ambient temperature or ambient specific enthalpy", enable=not 
            use_T_ambient));
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      " Ambient mass fractions m_i/m" annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nX > 0));
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100, 80; 100, -80], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-136, 144; 132, 82], string="%name")),
      Documentation(info="<html>
<p>
This element defines constant values for ambient pressure,
temperature and mass fractions. Note, that ambient temperature
and mass fractions have only an effect if the mass flow
is not, as usual, from the port in to the ambient,
but from the ambient in to the port.
</p>
</html>"));
    
  equation 
    Interfaces.checkAmbient(Medium.mediumName, Medium.incompressible,
      use_p_ambient, Medium.reducedX, Medium.nX, X_ambient);
    if use_p_ambient or Medium.incompressible then
      medium.p = p_ambient;
    else
      medium.d = d_ambient;
    end if;
    
    if use_T_ambient then
      medium.T = T_ambient;
    else
      medium.h = h_ambient;
    end if;
    
    medium.X = X_ambient;
  end FixedAmbient;
  
  model FixedMassFlowSource 
    "Ideal pump that produces a constant mass flow rate from a large reservoir" 
    
    import Modelica.SIunits.Conversions.*;
    
    extends Interfaces.PartialSource(medium(known_Th=if use_T_ambient then 
            Medium.Choices.Th.T_known else Medium.Choices.Th.h_known));
    
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate from an infinite reservoir in to the port";
    parameter Boolean use_T_ambient=true 
      "|Ambient temperature or ambient specific enthalpy| = true, if T_ambient is used, otherwise h_ambient"
      annotation (Evaluate=true);
    parameter Medium.Temperature T_ambient=from_degC(20) 
      "|Ambient temperature or ambient specific enthalpy| Ambient temperature, if use_T_ambient = true";
    parameter Medium.SpecificEnthalpy h_ambient=1.e4 
      "|Ambient temperature or ambient specific enthalpy| Ambient specific enthalpy, if use_T_ambient = false";
    parameter Medium.MassFraction X_ambient[Medium.nX](quantity=Medium.
          substanceNames) = zeros(Medium.nX) 
      "|Only for multi-substance flow| Ambient mass fractions m_i/m";
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(
        Rectangle(extent=[20, 60; 100, -60], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[38, 40; 100, -40], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
        Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
              fillColor=73)),
        Text(
          extent=[-54, 32; 16, -30],
          style(color=41, fillColor=41),
          string="m"),
        Text(extent=[-142, 142; 156, 88], string="%name"),
        Text(
          extent=[-124, -92; 148, -122],
          style(color=0),
          string="%m_flow"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1))),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Interfaces.checkAmbient(Medium.mediumName, Medium.incompressible, true,
      Medium.reducedX, Medium.nX, X_ambient);
    if use_T_ambient then
      medium.T = T_ambient;
    else
      medium.h = h_ambient;
    end if;
    medium.X = X_ambient;
    port.m_flow = -m_flow;
    
  end FixedMassFlowSource;
  
  model MassFlowSource 
    "Ideal pump that produces a mass flow rate from a large reservoir defined by input signal" 
    
    import Modelica.SIunits.Conversions.*;
    
    extends Interfaces.PartialSource(medium(known_Th=if use_T_ambient then 
            Medium.Choices.Th.T_known else Medium.Choices.Th.h_known));
    
    parameter Boolean use_T_ambient=true 
      "|Ambient temperature or ambient specific enthalpy| = true, if T_ambient is used, otherwise h_ambient"
      annotation (Evaluate=true);
    parameter Medium.Temperature T_ambient=from_degC(20) 
      "|Ambient temperature or ambient specific enthalpy| Ambient temperature, if use_T_ambient = true";
    parameter Medium.SpecificEnthalpy h_ambient=1.e4 
      "|Ambient temperature or ambient specific enthalpy| Ambient specific enthalpy, if use_T_ambient = false";
    parameter Medium.MassFraction X_ambient[Medium.nX](quantity=Medium.
          substanceNames) = zeros(Medium.nX) 
      "|Only for multi-substance flow| Ambient mass fractions m_i/m";
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(
        Rectangle(extent=[20, 60; 100, -60], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[38, 40; 100, -40], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
        Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
              fillColor=73)),
        Text(
          extent=[-54, 32; 16, -30],
          style(color=41, fillColor=41),
          string="m"),
        Text(extent=[-142, 142; 156, 88], string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1))),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
    Modelica.Blocks.Interfaces.RealInput m_flow(
                                            redeclare type SignalType = 
          Medium.MassFlowRate) 
      "Mass flow rate from an infinite reservoir in to the port as signal" 
      annotation (extent=[-140, -20; -100, 20]);
  equation 
    Interfaces.checkAmbient(Medium.mediumName, Medium.incompressible, true,
      Medium.reducedX, Medium.nX, X_ambient);
    if use_T_ambient then
      medium.T = T_ambient;
    else
      medium.h = h_ambient;
    end if;
    medium.X = X_ambient;
    port.m_flow = -m_flow;
  end MassFlowSource;
end Sources;
