import "C:/Dokumente und Einstellungen/proelss/Modelica/Subversion/Modelica_Fluid FreshFluid/Modelica_Fluid/package.mo";


package FluidConcept 
import SI = Modelica.SIunits;
  package Interfaces "Connectors" 
    connector AlphaFlowPortA 
      SI.Pressure p;
      flow SI.MassFlowRate m_flow;
      SI.SpecificEnthalpy h;
      flow SI.EnthalpyFlowRate H_flow;
      annotation (Icon(Ellipse(extent=[-60,80; 60,-40], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=69,
              rgbfillColor={0,128,255})), Text(
            extent=[-100,-20; 100,-100],
            string="%name",
            style(
              color=3,
              rgbcolor={0,0,255},
              pattern=0,
              fillColor=69,
              rgbfillColor={0,128,255},
              fillPattern=1))));
    end AlphaFlowPortA;
    
    connector AlphaFlowPortB 
      SI.Pressure p;
      flow SI.MassFlowRate m_flow;
      SI.SpecificEnthalpy h;
      flow SI.EnthalpyFlowRate H_flow;
      annotation (Icon(Ellipse(extent=[-60,80; 60,-40], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=69,
              rgbfillColor={0,128,255})),
          Ellipse(extent=[-40,60; 40,-20], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
                                          Text(
            extent=[-100,-20; 100,-100],
            string="%name",
            style(
              color=3,
              rgbcolor={0,0,255},
              pattern=0,
              fillColor=69,
              rgbfillColor={0,128,255},
              fillPattern=1))));
    end AlphaFlowPortB;
    
    connector BetaFlowPortA 
      extends PartialBetaFlowPort;
      input SI.SpecificEnthalpy h_ba "Property flow from b to a";
      output SI.SpecificEnthalpy h_ab "Property flow from a to b";
      annotation (Icon(Ellipse(extent=[-60,80; 60,-40], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=69,
              rgbfillColor={0,128,255})), Text(
            extent=[-100,-20; 100,-100],
            string="%name",
            style(
              color=3,
              rgbcolor={0,0,255},
              pattern=0,
              fillColor=69,
              rgbfillColor={0,128,255},
              fillPattern=1))));
    end BetaFlowPortA;
    
    connector BetaFlowPortB 
      extends PartialBetaFlowPort;
      input SI.SpecificEnthalpy h_ab "Property flow from a to b";
      output SI.SpecificEnthalpy h_ba "Property flow from b to a";
      annotation (Icon(
          Ellipse(extent=[-60,80; 60,-40], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=69,
              rgbfillColor={0,128,255})),
          Ellipse(extent=[-40,60; 40,-20], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
                                          Text(
            extent=[-100,-20; 100,-100],
            string="%name",
            style(
              color=3,
              rgbcolor={0,0,255},
              pattern=0,
              fillColor=69,
              rgbfillColor={0,128,255},
              fillPattern=1))));
    end BetaFlowPortB;
    
    partial connector PartialBetaFlowPort 
      SI.Pressure p;
      flow SI.MassFlowRate m_flow;
    end PartialBetaFlowPort;
  end Interfaces;
  annotation (uses(Modelica(version="2.2.2"), Modelica_Fluid(version=
            "1.0 Beta 3")), Icon, 
    Documentation(info="<html>
<p>This package provides basic component models to illustrate different approaches of 1D-fluid modelling which can then be compaired directly. The key issue which leads to the main differences of the approaches in this library is how to pass fluid properties, which are conveyed with the total massflow like specific enthalpy, mass fractions etc, across component boundaries. The different concepts also lead to variations in connector design.<p>
<p>In order to simplify the problem only specific enthalpy is conveyed with massflow, all other properties are neglected for demonstration puposes. The library is divided into subpackages which each contain a group of basic component models like control volumes, flow models or distributed flows, junctions and boundary conditions. In the Systems package similar system setups using different base concepts are collected.</p>
<p>
The following base concepts are featured in this collection:</p>
<b>Alpha</b>, current Modelica_Fluid approach:<br>
<ul>
<li>\"physical\" connector, contains specific enthalpy and enthalpy flow</li>
<li>\"semiLinear\"-operator handles changes in flow direction and zero mass flow, also allows for advanced optimization procedured on the tool side (common subexpression elimination)</li>
<li>port enthalpy depends on flow direction (may be undefined at zero flow) and connected components (either mixing enthalpy, if only components using semiLinear are connected or port volume enthalpy, if component with dynamic state at its port is attached) 
</ul>
<b>Beta</b>, ThermoPower-like approach:<br>
<ul>
<li>specific enthalpies on both sides of the connection set are included in the connector, no enthalpy flow, one-to-one restriction for connections in Modelica 3.0</li>
<li>change of flow direction is handled in balance equations</li>
<li>connector properties are always well defined</li>
</ul>
<b>Gamma</b>, Thermofluid-like approach:<br>
<ul>
<li>\"physical\" connector, contains specific enthalpy and enthalpy flow</li>
<li>change of flow direction is handled in balance equations</li>
<li>specific enthalpy is always passed to the \"left\", well defined port properties, efficient, but not very intuitive, restrictions on connections apply</li>
</ul>
The main problems all approaches have in common are:
<ul>
<li>port properties may not belong to the inside of the component and may easily be used in a wrong way</li>
<li>dynamic volumes are required to break large non-linear systems, especially in systems with non-linear flow models and lots of branches.</li>
<li>correct determination of the static head may require downstream properties inside the component boundary (not included yet)</li>
<li>Is a strict alternation of flow models and control volumes required? Is it possible to connect components in an arbitrary way? Keep pitfalls to a miminum.
</ul>


</html>", revisions="<html>
<ul>
<li><i>November 2007</i>
    by <a href=\"mailto:K.Proelss@tuhh.de\">Katrin Pr&ouml;l&szlig;</a>:<br>
    Created, based on discussions in the Fluid group</li>
<li><i>January 2008</i>
 by Katrin Pr&ouml;l&szlig;:<br>
    Added control volumes and flow models as well as a little bit of documentation</li>
</ul>
</html>"));
  package DistributedFlow 
    "Combines control volumes and flow models in a staggered grid" 
    partial model PartialDistributedFlow 
      replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium;
      replaceable package PressureDrop = WallFriction.Laminar extends 
        WallFriction.PartialWallFriction annotation(choicesAllMatching);
      parameter Integer n=10;
      parameter SI.Length length;
      parameter SI.Diameter diameter;
      parameter SI.Area Ah=Modelica.Constants.pi*diameter*length;
      parameter SI.Area Ac=Modelica.Constants.pi/4*diameter*diameter;
      parameter SI.Volume V=Modelica.Constants.pi/4*diameter*diameter*length;
      parameter SI.CoefficientOfHeatTransfer alpha0=300;
      Medium.BaseProperties[n] medium(each preferredMediumStates=true, p(start=p0),
                                    h(each start = h0), T(start=T0_init));
      SI.Mass[n] M;
      SI.InternalEnergy[n] U;
      SI.EnthalpyFlowRate[n + 1] H_flow;
      SI.MassFlowRate[n + 1] m_flow;
      SI.Pressure[n + 1] p;
      SI.Density[n] d_up;
      SI.DynamicViscosity[n] eta_up;
      annotation (Icon(
          Rectangle(extent=[-100,40; 100,-40], style(
              color=71,
              rgbcolor={85,170,255},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255})),
          Text(
            extent=[-80,-60; 80,-100],
            string="%name",
            style(
              color=3,
              rgbcolor={0,0,255},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}))),
                              Diagram);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] heatPort 
        annotation (extent=[-10,40; 10,60]);
    //Initialization
       parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.
            Init.NoInit "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean useInitialTemperature = true annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0_a annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0_b annotation(Dialog(tab = "Initialization"));
      parameter Medium.Temperature T0 annotation(Dialog(tab = "Initialization",enable=useInitialTemperature));
      parameter Medium.SpecificEnthalpy h0 = Medium.specificEnthalpy(Medium.setState_pTX(0.5*(p0_a+p0_b),T0)) 
                                                                              annotation(Dialog(tab = "Initialization",enable=not useInitialTemperature));
    protected 
      parameter Medium.AbsolutePressure[n] p0=if n>1 then linspace(p0_a,p0_b,n) else fill(0.5*(p0_a+p0_b),n);
      parameter Medium.Temperature[n] T0_init=if useInitialTemperature then fill(T0,n) else {Medium.temperature(Medium.setState_phX(p0[i],h0)) for i in 1:n};
    initial equation 
      if not Medium.singleState then
        if initType==Modelica_Fluid.Types.Init.SteadyState then
          der(medium.p) = zeros(n);
        elseif initType==Modelica_Fluid.Types.Init.InitialValues then
          medium.p = p0;
        end if;
      end if;
      if initType==Modelica_Fluid.Types.Init.SteadyState then
        medium.h = zeros(n);
      elseif initType==Modelica_Fluid.Types.Init.InitialValues then
        medium.h = fill(h0,n);
      end if;
    equation 
    //Properties depending on flow direction and on dynamic states
      for i in 2:n loop
        H_flow[i] = if m_flow[i] >= 0 then m_flow[i]*medium[i - 1].h else m_flow[i]
          *medium[i].h;
      d_up[i - 1] = if m_flow[i] >= 0 then medium[i - 1].d else medium[i].d;
       //d_up[i-1] = medium[i-1].d;
       eta_up[i - 1] = if m_flow[i] >= 0 then Medium.dynamicViscosity(medium[i - 1].state) else 
                Medium.dynamicViscosity(medium[i].state);
     //eta_up[i-1]=Medium.dynamicViscosity(medium[i - 1].state);
      end for;
      
      for i in 1:n loop
        M[i] = medium[i].d*V/n;
        U[i] = medium[i].u*M[i];
    //Energy balance
        der(U[i]) = H_flow[i] - H_flow[i + 1] + heatPort[i].Q_flow;
    //Mass balance
        der(M[i]) = m_flow[i] - m_flow[i + 1];
      end for;
      
    //Pressure Drop: mflow=f(dp)
      for i in 1:n loop
        medium[i].p = p[i];
        m_flow[i+1]=PressureDrop.massFlowRate_dp(
        p[i] - p[i + 1],
          d_up[i],
          eta_up[i],
          length/n,
          diameter,
          roughness=2.5e-5,
          dp_small=1);
      end for;
      
    //HeatTransfer
    for i in 1:n loop
      heatPort[i].Q_flow=alpha0*Ah/n*(heatPort[i].T-medium[i].T);
    end for;
      
    end PartialDistributedFlow;
    
    model AlphaDistributedFlow 
      extends FluidConcept.DistributedFlow.PartialDistributedFlow;
      FluidConcept.Interfaces.AlphaFlowPortA port_a 
                                      annotation (extent=[-124,-18; -94,12]);
      FluidConcept.Interfaces.AlphaFlowPortB port_b 
                                      annotation (extent=[94,-18; 124,12]);
      annotation (Icon(Text(
            extent=[-68,12; 60,-12],
            string="alpha",
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}))));
    equation 
      
    //Port variables
      port_a.p=medium[1].p;
      port_b.p=p[n+1];
      port_a.m_flow=m_flow[1];
      port_b.m_flow=-m_flow[n+1];
      port_a.H_flow=H_flow[1];
      port_b.H_flow=-H_flow[n+1];
      
    //Interface properties
      port_a.H_flow = semiLinear(
        port_a.m_flow,
        port_a.h,
        medium[1].h);
      port_b.H_flow = semiLinear(
        port_b.m_flow,
        port_b.h,
        medium[n].h);
       d_up[n] = if m_flow[n + 1] >= 0 then medium[n].d else Medium.density_phX(
        port_b.p, port_b.h);
       eta_up[n] = if m_flow[n + 1] >= 0 then Medium.dynamicViscosity(medium[n].state) else 
              Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h));
      // d_up[n] =  medium[n].d;
      // eta_up[n] = Medium.dynamicViscosity(medium[n].state);
    end AlphaDistributedFlow;
    
    model BetaDistributedFlow 
      extends FluidConcept.DistributedFlow.PartialDistributedFlow;
      FluidConcept.Interfaces.BetaFlowPortA port_a 
                                      annotation (extent=[-124,-18; -94,12]);
      FluidConcept.Interfaces.BetaFlowPortB port_b 
                                      annotation (extent=[94,-18; 124,12]);
      annotation (Icon(Text(
            extent=[-62,12; 66,-12],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="beta")));
    equation 
      
    //Port variables
      port_a.p=medium[1].p;
      port_b.p=p[n+1];
      port_a.m_flow=m_flow[1];
      port_b.m_flow=-m_flow[n+1];
      
      port_a.h_ab=medium[1].h;
      port_b.h_ba=medium[n].h;
      
      H_flow[1]=smooth(0,if m_flow[1]>=0 then m_flow[1]*port_a.h_ba else port_a.m_flow*port_a.h_ab);
      H_flow[n+1]=smooth(0,if m_flow[n+1]<=0 then m_flow[n+1]*port_b.h_ab else m_flow[n+1]*port_b.h_ba);
      d_up[n] = if m_flow[n + 1] >= 0 then medium[n].d else Medium.density_phX(
        port_b.p, port_b.h_ab);
      eta_up[n] = if m_flow[n + 1] >= 0 then Medium.dynamicViscosity(medium[n].state) else 
              Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h_ab));
      
    end BetaDistributedFlow;
    
    model GammaDistributedFlow 
      extends FluidConcept.DistributedFlow.PartialDistributedFlow;
      FluidConcept.Interfaces.AlphaFlowPortA port_a 
                                      annotation (extent=[-124,-18; -94,12]);
      Interfaces.AlphaFlowPortB port_b 
                                      annotation (extent=[94,-18; 124,12]);
      annotation (Icon(Text(
            extent=[-64,12; 64,-12],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="gamma")));
    equation 
      
    //Port variables
      port_a.p=medium[1].p;
      port_b.p=p[n+1];
      port_a.m_flow=m_flow[1];
      port_b.m_flow=-m_flow[n+1];
      port_a.H_flow=H_flow[1];
      port_b.H_flow=-H_flow[n+1];
      port_a.h=medium[1].h;
      port_b.H_flow=noEvent(if port_b.m_flow >= 0 then port_b.m_flow*port_b.h else port_b.m_flow*medium[n].h);
      
    //upstream properties at port_b
      d_up[n] = if m_flow[n + 1] >= 0 then medium[n].d else Medium.density_phX(
        port_b.p, port_b.h);
      eta_up[n] = if m_flow[n + 1] >= 0 then Medium.dynamicViscosity(medium[n].state) else 
              Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h));
      
    end GammaDistributedFlow;
    
    model DeltaDistributedFlow 
      replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium;
      replaceable package PressureDrop = WallFriction.Laminar extends 
        WallFriction.PartialWallFriction 
                                       annotation(choicesAllMatching);
      parameter Integer n=10;
      parameter SI.Length length;
      parameter SI.Diameter diameter;
      parameter SI.Area Ah=Modelica.Constants.pi*diameter*length;
      parameter SI.Area Ac=Modelica.Constants.pi/4*diameter*diameter;
      parameter SI.Volume V=Modelica.Constants.pi/4*diameter*diameter*length;
      parameter SI.CoefficientOfHeatTransfer alpha0=300;
      Medium.BaseProperties[n] medium(
        each preferredMediumStates=true,
        p(start=p0),
        h(each start=h0),
        T(start=T0_init));
      SI.Mass[n] M;
      SI.InternalEnergy[n] U;
      SI.EnthalpyFlowRate[n + 1] H_flow;
      SI.MassFlowRate[n + 1] m_flow;
      //SI.Pressure[n + 2] p;
      SI.Density[n + 1] d_up(each start=Medium.density(Medium.setState_phX(0.5*(p0_a+p0_b),h0)));
      SI.DynamicViscosity[n + 1] eta_up(each start=Medium.dynamicViscosity(Medium.setState_phX(0.5*(p0_a+p0_b),h0)));
      annotation (Icon(
          Rectangle(extent=[-100,40; 100,-40], style(
              color=71,
              rgbcolor={85,170,255},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255})),
          Text(
            extent=[-80,-60; 80,-100],
            string="%name",
            style(
              color=3,
              rgbcolor={0,0,255},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255})),
          Text(
            extent=[-64,12; 64,-12],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="delta")), Diagram);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] heatPort 
        annotation (extent=[-10,40; 10,60]);
    //Initialization
      parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.Init.NoInit 
        "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean useInitialTemperature=true   annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0_a annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0_b annotation(Dialog(tab = "Initialization"));
      parameter Medium.Temperature T0 annotation(Dialog(tab = "Initialization",enable=useInitialTemperature));
      parameter Medium.SpecificEnthalpy h0=Medium.specificEnthalpy(
          Medium.setState_pTX(0.5*(p0_a + p0_b), T0))                         annotation(Dialog(tab = "Initialization",enable=not useInitialTemperature));
      
      Interfaces.AlphaFlowPortA port_a(p(start=p0_a),h(start=h0)) 
                                      annotation (extent=[-116,-16; -86,14]);
      Interfaces.AlphaFlowPortA port_b(p(start=p0_b),h(start=h0)) 
                                      annotation (extent=[94,-18; 124,12]);
      annotation (Icon);
    protected 
      parameter Medium.AbsolutePressure[n] p0={(p0_b - p0_a)/(2*n)*(2*(i - 1) + 1)
           + p0_a for i in 1:n};
      parameter Medium.Temperature[n] T0_init=if useInitialTemperature then fill(T0,
          n) else {Medium.temperature(Medium.setState_phX(p0[i], h0)) for i in 1:n};
    initial equation 
      if not Medium.singleState then
        if initType == Modelica_Fluid.Types.Init.SteadyState then
          der(medium.p) = zeros(n);
        elseif initType == Modelica_Fluid.Types.Init.InitialValues then
          medium.p = p0;
        end if;
      end if;
      if initType == Modelica_Fluid.Types.Init.SteadyState then
        medium.h = zeros(n);
      elseif initType == Modelica_Fluid.Types.Init.InitialValues then
        medium.h = fill(h0, n);
      end if;
      
    equation 
      
    //Properties depending on flow direction and on dynamic states
      for i in 2:n loop
        H_flow[i] = if m_flow[i] >= 0 then m_flow[i]*medium[i - 1].h else m_flow[i]
          *medium[i].h;
      d_up[i] = if m_flow[i] >= 0 then medium[i - 1].d else medium[i].d;
      //  d_up[i] = medium[i - 1].d;
        eta_up[i] = if m_flow[i] >= 0 then Medium.dynamicViscosity(medium[i - 1].state) else 
                 Medium.dynamicViscosity(medium[i].state);
    //    eta_up[i] = Medium.dynamicViscosity(medium[i - 1].state);
      end for;
      
      for i in 1:n loop
        M[i] = medium[i].d*V/n;
        U[i] = medium[i].u*M[i];
    //Energy balance
        der(U[i]) = H_flow[i] - H_flow[i + 1] + heatPort[i].Q_flow;
    //Mass balance
        der(M[i]) = m_flow[i] - m_flow[i + 1];
      end for;
      
    //Pressure Drop: mflow=f(dp)
      
      for i in 1:n - 1 loop
        m_flow[i + 1] = PressureDrop.massFlowRate_dp(
          medium[i].p - medium[i + 1].p,
          d_up[i + 1],
          eta_up[i + 1],
          length/n,
          diameter,
          roughness=2.5e-5,
          dp_small=1);
      end for;
      
       m_flow[1] = PressureDrop.massFlowRate_dp(
        (port_a.p - medium[1].p),
        d_up[1],
        eta_up[1],
        length/(2*n),
        diameter,
        roughness=2.5e-5,
        dp_small=1);
      m_flow[n + 1] = PressureDrop.massFlowRate_dp(
        (medium[n].p - port_b.p),
        d_up[n + 1],
         eta_up[n + 1],
         length/(2*n),
        diameter,
        roughness=2.5e-5,
       dp_small=1);
      
    //HeatTransfer
      for i in 1:n loop
        heatPort[i].Q_flow = alpha0*Ah/n*(heatPort[i].T - medium[i].T);
      end for;
      
    //Port variables
      port_a.m_flow = m_flow[1];
      port_b.m_flow = -m_flow[n + 1];
      port_a.H_flow = H_flow[1];
      port_b.H_flow = -H_flow[n + 1];
      
    //Interface properties
      port_a.H_flow = semiLinear(
        port_a.m_flow,
        port_a.h,
        medium[1].h);
      port_b.H_flow = semiLinear(
        port_b.m_flow,
        port_b.h,
        medium[n].h);
       d_up[n+1] = if m_flow[n + 1] >= 0 then medium[n].d else Medium.density_phX(
        port_b.p, port_b.h);
       d_up[1] = if m_flow[1] >= 0 then Medium.density_phX(
         port_a.p, port_a.h) else medium[1].d;
       eta_up[1] = if m_flow[1] >= 0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p, port_a.h)) else 
                          Medium.dynamicViscosity(medium[1].state);
       eta_up[n+1] = if m_flow[n + 1] >= 0 then Medium.dynamicViscosity(medium[n].state) else 
              Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h));
      // d_up[1] = medium[1].d;
      // eta_up[1] = Medium.dynamicViscosity(medium[1].state);
      // d_up[n + 1] = medium[n].d;
      // eta_up[n + 1] = Medium.dynamicViscosity(medium[n].state);
    end DeltaDistributedFlow;
  end DistributedFlow;
  
  package ControlVolumes "Dynamic control volumes" 
    partial model PartialCV 
    replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium  annotation(choicesAllMatching);
    parameter SI.Length length;
    parameter SI.Diameter diameter;
    parameter SI.Area Ah=Modelica.Constants.pi*diameter*length;
    parameter SI.Area Ac=Modelica.Constants.pi/4*diameter*diameter;
    parameter SI.Volume V=Modelica.Constants.pi/4*diameter*diameter*length;
    parameter SI.CoefficientOfHeatTransfer alpha0=300;
    Medium.BaseProperties medium(preferredMediumStates=true, p(start=p0),
                                  h(start = h0), T(start=T0_init));
    SI.Mass M;
    SI.InternalEnergy U;
    SI.EnthalpyFlowRate H_flow;
    SI.MassFlowRate m_flow;
      
    annotation (Icon(
        Text(
          extent=[-80,-60; 80,-100],
          string="%name",
          style(
            color=3,
            rgbcolor={0,0,255},
            gradient=2,
            fillColor=71,
            rgbfillColor={85,170,255})), Ellipse(extent=[-60,60; 60,-60], style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}))),
                            Diagram);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (extent=[-10,60; 10,80]);
      //Initialization
     parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Boolean useInitialTemperature = true annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p0 annotation(Dialog(tab = "Initialization"));
    parameter Medium.Temperature T0 annotation(Dialog(tab = "Initialization",enable=useInitialTemperature));
    parameter Medium.SpecificEnthalpy h0 = Medium.specificEnthalpy(Medium.setState_pTX(p0,T0)) 
                                                                            annotation(Dialog(tab = "Initialization",enable=not useInitialTemperature));
    protected 
     parameter Medium.Temperature T0_init=if useInitialTemperature then T0 else Medium.temperature(Medium.setState_phX(p0,h0));
    initial equation 
    if not Medium.singleState then
      if initType==Modelica_Fluid.Types.Init.SteadyState then
        der(medium.p) = 0;
      elseif initType==Modelica_Fluid.Types.Init.InitialValues then
        medium.p = p0;
      end if;
    end if;
    if initType==Modelica_Fluid.Types.Init.SteadyState then
      medium.h = 0;
    elseif initType==Modelica_Fluid.Types.Init.InitialValues then
      medium.h =h0;
    end if;
    equation 
      M = medium.d*V;
      U = medium.u*M;
        //Energy balance
      der(U) = H_flow + heatPort.Q_flow;
        //Mass balance
      der(M) = m_flow;
      
      //HeatTransfer
      heatPort.Q_flow=alpha0*Ah*(heatPort.T-medium.T);
      
    end PartialCV;
    
    model AlphaCVPortVolume 
      extends FluidConcept.ControlVolumes.PartialCV;
      Interfaces.AlphaFlowPortA port annotation (extent=[-14,-18; 16,12]);
      annotation (Icon(Text(
            extent=[-64,34; 64,10],
            string="alpha",
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}))));
    equation 
      H_flow=port.H_flow;
      m_flow=port.m_flow;
      port.h=medium.h;
      port.p=medium.p;
    end AlphaCVPortVolume;
    
    model AlphaCV 
      extends FluidConcept.ControlVolumes.PartialCV;
      Interfaces.AlphaFlowPortB port_a 
                                     annotation (extent=[-84,-18; -54,12]);
      annotation (Icon(Text(
            extent=[-62,14; 66,-10],
            string="alpha",
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}))));
      Interfaces.AlphaFlowPortB port_b annotation (extent=[54,-18; 84,12]);
    equation 
      port_a.H_flow=semiLinear(port_a.m_flow,port_a.h,medium.h);
      port_b.H_flow=semiLinear(port_b.m_flow,port_b.h,medium.h);
      H_flow=port_a.H_flow+port_b.H_flow;
      m_flow=port_a.m_flow+port_b.m_flow;
      port_a.p=medium.p;
      port_b.p=medium.p;
    end AlphaCV;
    
    model BetaCV 
      extends FluidConcept.ControlVolumes.PartialCV;
      parameter Integer na=1 "number of a ports";
      parameter Integer nb=1 "number of b ports";
      Interfaces.BetaFlowPortA[na] port_a 
                                      annotation (extent=[-84,-18; -54,12]);
      annotation (Icon(Text(
            extent=[-64,14; 64,-10],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="beta")));
      Interfaces.BetaFlowPortB[nb] port_b 
                                      annotation (extent=[54,-18; 84,12]);
    equation 
      m_flow=sum(port_a.m_flow)+sum(port_b.m_flow);
      H_flow=sum({if port_a[i].m_flow>=0 then port_a[i].m_flow*port_a[i].h_ba else port_a[i].m_flow*medium.h for i in 1:na})
            +sum({if port_b[i].m_flow>=0 then port_b[i].m_flow*port_b[i].h_ab else port_b[i].m_flow*medium.h for i in 1:nb});
      port_a.p=fill(medium.p,na);
      port_b.p=fill(medium.p,nb);
      port_a.h_ab=fill(medium.h,na);
      port_b.h_ba=fill(medium.h,nb);
    end BetaCV;
    
    model GammaCV 
      extends FluidConcept.ControlVolumes.PartialCV;
      Interfaces.AlphaFlowPortA port   annotation (extent=[-14,-18; 16,12]);
    equation 
      port.h=medium.h;
      H_flow=port.H_flow;
      m_flow=port.m_flow;
      port.p=medium.p;
      
      annotation (Icon(Text(
            extent=[-62,38; 66,14],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="gamma")), Diagram);
    end GammaCV;
  end ControlVolumes;
  
  package FlowModels "Algebraic flow models" 
    partial model PartialFlowModel 
      replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
      replaceable package PressureDrop = WallFriction.Laminar extends 
        WallFriction.PartialWallFriction   annotation(choicesAllMatching);
      parameter SI.Length length;
      parameter SI.Diameter diameter;
      parameter Medium.AbsolutePressure p0_a annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0_b annotation(Dialog(tab = "Initialization"));
      parameter SI.MassFlowRate m_flow0 annotation(Dialog(tab = "Initialization"));
      SI.Pressure[2] p(start={p0_a,p0_b});
      SI.MassFlowRate[2] m_flow(each start=m_flow0) "port massflows";
      SI.Density d_up;
      SI.DynamicViscosity eta_up;
    equation 
    //Mass balance
      0=m_flow[1]-m_flow[2];
    //Pressure Drop: mflow=f(dp)
          m_flow[2]=PressureDrop.massFlowRate_dp(
          p[1] - p[2],
          d_up,
          eta_up,
          length,
          diameter,
          roughness=2.5e-5,
          dp_small=1);
      annotation (Icon(Rectangle(extent=[-60,60; 60,-40], style(
              color=7,
              rgbcolor={255,255,255},
              thickness=4,
              fillColor=7,
              rgbfillColor={255,255,255})), Line(points=[-60,-40; -60,60; 60,-40;
                60,60], style(
              color=3,
              rgbcolor={0,0,255},
              thickness=4))));
    end PartialFlowModel;
    
    model AlphaFM 
      extends PartialFlowModel;
      annotation (Icon(
        Text(
          extent=[-74,-40; 86,-80],
          string="%name",
          style(
            color=3,
            rgbcolor={0,0,255},
            gradient=2,
            fillColor=71,
            rgbfillColor={85,170,255})),
                       Text(
            extent=[-48,64; 80,40],
            string="alpha",
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}))));
      Interfaces.AlphaFlowPortB port_a annotation (extent=[-84,-8; -54,22]);
      Interfaces.AlphaFlowPortB port_b annotation (extent=[54,-8; 84,22]);
    equation 
      port_a.p=p[1];
      port_b.p=p[2];
      port_a.m_flow=m_flow[1];
      port_b.m_flow=-m_flow[2];
      port_a.H_flow=semiLinear(port_a.m_flow,port_a.h,port_b.h);
      0=port_a.H_flow+port_b.H_flow;
      d_up=if port_a.m_flow>=0 then Medium.density(Medium.setState_phX(port_a.p,port_a.h)) else 
        Medium.density(Medium.setState_phX(port_b.p,port_b.h));
      eta_up=if port_a.m_flow>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p,port_a.h)) else 
        Medium.dynamicViscosity(Medium.setState_phX(port_b.p,port_b.h));
    end AlphaFM;
    
    model BetaFM "isenthalpic flow model" 
      extends PartialFlowModel;
      annotation (Icon(
        Text(
          extent=[-74,-40; 86,-80],
          string="%name",
          style(
            color=3,
            rgbcolor={0,0,255},
            gradient=2,
            fillColor=71,
            rgbfillColor={85,170,255})),
                       Text(
            extent=[-46,52; 82,28],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="beta")));
      Interfaces.BetaFlowPortA port_a annotation (extent=[-84,-8; -54,22]);
      Interfaces.BetaFlowPortB port_b annotation (extent=[54,-8; 84,22]);
    equation 
      port_a.p=p[1];
      port_b.p=p[2];
      port_a.m_flow=m_flow[1];
      port_b.m_flow=-m_flow[2];
      port_a.h_ab=port_b.h_ab;
      port_a.h_ba=port_b.h_ba;
      d_up=if port_a.m_flow>=0 then Medium.density(Medium.setState_phX(port_a.p,port_a.h_ba)) else 
        Medium.density(Medium.setState_phX(port_b.p,port_b.h_ab));
      eta_up=if port_a.m_flow>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p,port_a.h_ba)) else 
        Medium.dynamicViscosity(Medium.setState_phX(port_b.p,port_b.h_ab));
    end BetaFM;
    
    model GammaFM 
      extends PartialFlowModel;
      annotation (Icon(
        Text(
          extent=[-74,-40; 86,-80],
          string="%name",
          style(
            color=3,
            rgbcolor={0,0,255},
            gradient=2,
            fillColor=71,
            rgbfillColor={85,170,255})),
                       Text(
            extent=[-46,78; 82,54],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="gamma")));
      Interfaces.AlphaFlowPortB port_a annotation (extent=[-84,-8; -54,22]);
      Interfaces.AlphaFlowPortB port_b annotation (extent=[54,-8; 84,22]);
    equation 
      port_a.p=p[1];
      port_b.p=p[2];
      port_a.m_flow=m_flow[1];
      port_b.m_flow=-m_flow[2];
      port_a.H_flow=if port_a.m_flow>=0 then port_a.m_flow*port_a.h else port_a.m_flow*port_b.h;
      port_a.H_flow=-port_b.H_flow;
      d_up=if port_a.m_flow>=0 then Medium.density(Medium.setState_phX(port_a.p,port_a.h)) else 
        Medium.density(Medium.setState_phX(port_b.p,port_b.h));
      eta_up=if port_a.m_flow>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p,port_a.h)) else 
        Medium.dynamicViscosity(Medium.setState_phX(port_b.p,port_b.h));
    end GammaFM;
  end FlowModels;
  
  package WallFriction 
    "Pressure drop package from Modelica_Fluid, adaptedDifferent variants for pressure drops due to pipe wall friction" 
    partial package PartialWallFriction 
      "Partial wall friction characteristic (base package of all wall friction characteristics)" 
      
      annotation (Documentation(info="<html>
 
</html>"));
      
    // Constants to be set in subpackages
      constant Boolean use_eta = true 
        "= true, if eta_a/eta_b are used in function, otherwise value is not used";
      constant Boolean use_roughness = true 
        "= true, if roughness is used in function, otherwise value is not used";
      constant Boolean use_dp_small = true 
        "= true, if dp_small is used in function, otherwise value is not used";
      constant Boolean use_m_flow_small = true 
        "= true, if m_flow_small is used in function, otherwise value is not used";
      constant Boolean dp_is_zero = false 
        "= true, if no wall friction is present, i.e., dp = 0 (function massFlowRate_dp() cannot be used)";
      
    // pressure loss characteristic functions
      replaceable partial function massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
                extends Modelica.Icons.Function;
        
        input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        input SI.Density d_up "Upstream density";
       // input SI.Density d_down "Downstream density";
        input SI.DynamicViscosity eta_up 
          "Upstream viscosity (dummy if use_eta = false)";
       // input SI.DynamicViscosity eta_down 
       //   "Downstream viscosity (dummy if use_eta = false)";
        input SI.Length length "Length of pipe";
        input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
        input SI.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
        input SI.AbsolutePressure dp_small = 1 
          "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";
        
        output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
      annotation (Documentation(info="<html>
 
</html>"));
      end massFlowRate_dp;
      
      replaceable partial function pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
                extends Modelica.Icons.Function;
        
        input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        input SI.Density d_up "Upstream density";
        input SI.Density d_down "Downstream density";
        input SI.DynamicViscosity eta_up 
          "Upstream viscosity (dummy if use_eta = false)";
        input SI.DynamicViscosity eta_down 
          "Downstream viscosity (dummy if use_eta = false)";
        // input SI.Density d_a "Density at port_a";
        // input SI.Density d_b "Density at port_b";
        // input SI.DynamicViscosity eta_a 
        //   "Dynamic viscosity at port_a (dummy if use_eta = false)";
        // input SI.DynamicViscosity eta_b 
        //   "Dynamic viscosity at port_b (dummy if use_eta = false)";
        input SI.Length length "Length of pipe";
        input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
        input SI.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
        input SI.MassFlowRate m_flow_small = 0.01 
          "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
        output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        
      annotation (Documentation(info="<html>
 
</html>"));
      end pressureLoss_m_flow;
      
    end PartialWallFriction;
    
    annotation (Documentation(info="<html>
<p>
This package provides functions to compute
pressure losses due to <b>wall friction</b> in a pipe.
Every correlation is defined by a package that is derived
by inheritance from the package WallFriction.PartialWallFriction.
The details of the underlying pipe wall friction model are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
Basically, different variants of the equation
</p>
 
<pre>
   dp = &lambda;(Re,<font face=\"Symbol\">D</font>)*(L/D)*&rho;*v*|v|/2
</pre>
 
<p>
are used, where the friction loss factor &lambda; is shown
in the next figure:
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));
    package NoFriction "No pipe wall friction" 
      
      annotation (Documentation(info="<html>
<p>
This component sets the pressure loss due to wall friction 
to zero, i.e., it allows to switch off pipe wall friction.
</p>
</html>"));
      
      extends PartialWallFriction(
                final use_eta = false,
                final use_roughness = false,
                final use_dp_small = false,
                final use_m_flow_small = false,
                final dp_is_zero = true);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        
        annotation (Documentation(info="<html>
 
</html>"));
      algorithm 
        assert(false, "function massFlowRate_dp (option: from_dp=true)
cannot be used for WallFriction.NoFriction. Use instead
function pressureLoss_m_flow (option: from_dp=false)");
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        
        annotation (Documentation(info="<html>
 
</html>"));
      algorithm 
        dp := 0;
      end pressureLoss_m_flow;
    end NoFriction;
    
    package Laminar 
      "Pipe wall friction in the laminar regime (linear correlation)" 
      
      annotation (Documentation(info="<html>
<p>
This component defines only the laminar region of wall friction:
dp = k*m_flow, where \"k\" depends on density and dynamic viscosity.
The roughness of the wall does not have an influence on the laminar
flow and therefore argument roughness is ignored.
Since this is a linear relationship, the occuring systems of equations
are usually much simpler (e.g. either linear instead of non-linear).
By using nominal values for density and dynamic viscosity, the 
systems of equations can still further be reduced. 
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the \"light blue curve\" called
<b>Hagen-Poiseuille</b>.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));
      
      extends PartialWallFriction(
                final use_eta = true,
                final use_roughness = false,
                final use_dp_small = false,
                final use_m_flow_small = false);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        
        annotation (Documentation(info="<html>
 
</html>"));
      algorithm 
       // m_flow :=dp*Modelica.Constants.pi*diameter^4*(d_a + d_b)/(128*length*(eta_a + eta_b));
       m_flow :=dp*Modelica.Constants.pi*diameter^4*d_up/(128*length*eta_up);
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        
        annotation (Documentation(info="<html>
 
</html>"));
      algorithm 
        //dp := m_flow*128*length*(eta_a+eta_b)/(Modelica.Constants.pi*diameter^4*(d_a + d_b));
        dp := m_flow*128*length*eta_up/(Modelica.Constants.pi*diameter^4*d_up);
      end pressureLoss_m_flow;
    end Laminar;
    
    package QuadraticTurbulent 
      "Pipe wall friction in the quadratic turbulent regime (simple characteristic, eta not used)" 
      
      annotation (Documentation(info="<html>
<p>
This component defines only the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is no longer a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the asymptotic behaviour for large
Reynolds numbers, i.e., the values at the right ordinate where
&lambda; is constant.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));
      
      extends PartialWallFriction(
                final use_eta = false,
                final use_roughness = true,
                final use_dp_small = true,
                final use_m_flow_small = true);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import Modelica.Math;
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi = Modelica.Constants.pi;
        Real zeta;
        Real k_inv;
      algorithm 
        /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
        assert(roughness > 1.e-10,
               "roughness > 0 required for quadratic turbulent wall friction characteristic");
        zeta  := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
        k_inv := (pi*diameter*diameter)^2/(8*zeta);
        m_flow := Modelica_Fluid.Utilities.regRoot2(dp, dp_small, d_up*k_inv, d_up*k_inv);
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        import Modelica.Math;
        
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi = Modelica.Constants.pi;
        Real zeta;
        Real k;
      algorithm 
        /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
        assert(roughness > 1.e-10,
               "roughness > 0 required for quadratic turbulent wall friction characteristic");
        zeta := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
        k    := 8*zeta/(pi*diameter*diameter)^2;
        dp   := Modelica_Fluid.Utilities.regSquare2(m_flow, m_flow_small, k/d_up, k/d_up);
      end pressureLoss_m_flow;
    end QuadraticTurbulent;
    
    package LaminarAndQuadraticTurbulent 
      "Pipe wall friction in the laminar and quadratic turbulent regime (simple characteristic)" 
      
      annotation (Documentation(info="<html>
<p>
This component defines the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is no longer a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
At Re=4000, a polynomial is constructed that approaches
the constant &lambda; (for large Reynolds-numbers) at Re=4000
smoothly and has a derivative at zero mass flow rate that is
identical to laminar wall friction.
</p>
</html>"));
      
      extends PartialWallFriction(
                final use_eta = true,
                final use_roughness = true,
                final use_dp_small = false,
                final use_m_flow_small = false);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import Modelica.Math;
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
        constant Real Re_turbulent = 4000 "Start of turbulent regime";
        Real zeta;
        Real k0;
        Real k_inv;
        Real yd0 "Derivative of m_flow=m_flow(dp) at zero";
        SI.AbsolutePressure dp_turbulent;
      algorithm 
      /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If data.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
        assert(roughness > 1.e-10,
               "roughness > 0 required for quadratic turbulent wall friction characteristic");
        zeta   := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
        k0     := 128*length/(pi*diameter^4);
        k_inv  := (pi*diameter*diameter)^2/(8*zeta);
        yd0    := d_up/(k0*eta_up);
        dp_turbulent := (eta_up*diameter*pi/4)^2*Re_turbulent^2/(k_inv*d_up);
        m_flow := Modelica_Fluid.Utilities.regRoot2(dp, dp_turbulent, d_up*k_inv, d_up*k_inv,
                                                    use_yd0=true, yd0=yd0);
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        import Modelica.Math;
        
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
        constant Real Re_turbulent = 4000 "Start of turbulent regime";
        Real zeta;
        Real k0;
        Real k;
        Real yd0 "Derivative of dp = f(m_flow) at zero";
        SI.MassFlowRate m_flow_turbulent 
          "The turbulent region is: |m_flow| >= m_flow_turbulent";
        
      algorithm 
      /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
*/
        assert(roughness > 1.e-10,
               "roughness > 0 required for quadratic turbulent wall friction characteristic");
        zeta := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
        k0   := 128*length/(pi*diameter^4);
        k    := 8*zeta/(pi*diameter*diameter)^2;
        yd0  := k0*eta_up/d_up;
        m_flow_turbulent :=(pi/4)*diameter*eta_up*Re_turbulent;
        dp :=Modelica_Fluid.Utilities.regSquare2(m_flow, m_flow_turbulent, k/d_up, k/d_up,
                                                 use_yd0=true, yd0=yd0);
      end pressureLoss_m_flow;
    end LaminarAndQuadraticTurbulent;
    
    package Detailed 
      "Pipe wall friction in the whole regime (detailed characteristic)" 
      
      annotation (Documentation(info="<html>
<p>
This component defines the complete regime of wall friction.
The details are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
The functional relationship of the friction loss factor &lambda; is
displayed in the next figure. Function massFlowRate_dp() defines the \"red curve\"
(\"Swamee and Jain\"), where as function pressureLoss_m_flow() defines the
\"blue curve\" (\"Colebrook-White\"). The two functions are inverses from 
each other and give slightly different results in the transition region
between Re = 1500 .. 4000, in order to get explicit equations without
solving a non-linear equation.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
</html>"));
      
      extends PartialWallFriction(
                final use_eta = true,
                final use_roughness = false,
                final use_dp_small = false,
                final use_m_flow_small = false);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import Modelica.Math;
                annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi = Modelica.Constants.pi;
        Real Delta = roughness/diameter "Relative roughness";
        SI.ReynoldsNumber Re1 = (745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97 
          "Re leaving laminar curve";
        SI.ReynoldsNumber Re2 = 4000 "Re entering turbulent curve";
        SI.DynamicViscosity eta "Upstream viscosity";
        SI.Density d "Upstream density";
        SI.ReynoldsNumber Re "Reynolds number";
        Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
        
        function interpolateInRegion2 
           input Real Re_turbulent;
           input SI.ReynoldsNumber Re1;
           input SI.ReynoldsNumber Re2;
           input Real Delta;
           input Real lambda2;
           output SI.ReynoldsNumber Re;
           annotation(smoothOrder=1);
          // point lg(lambda2(Re1)) with derivative at lg(Re1)
        protected 
          Real x1=Math.log10(64*Re1);
          Real y1=Math.log10(Re1);
          Real yd1=1;
          
          // Point lg(lambda2(Re2)) with derivative at lg(Re2)
          Real aux1=(0.5/Math.log(10))*5.74*0.9;
          Real aux2=Delta/3.7 + 5.74/Re2^0.9;
          Real aux3=Math.log10(aux2);
          Real L2=0.25*(Re2/aux3)^2;
          Real aux4=2.51/sqrt(L2) + 0.27*Delta;
          Real aux5=-2*sqrt(L2)*Math.log10(aux4);
          Real x2=Math.log10(L2);
          Real y2=Math.log10(aux5);
          Real yd2=0.5 + (2.51/Math.log(10))/(aux5*aux4);
          
          // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
          Real diff_x=x2 - x1;
          Real m=(y2 - y1)/diff_x;
          Real c2=(3*m - 2*yd1 - yd2)/diff_x;
          Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
          Real lambda2_1=64*Re1;
          Real dx;
        algorithm 
           dx := Math.log10(lambda2/lambda2_1);
           Re := Re1*(lambda2/lambda2_1)^(1 + dx*(c2 + dx*c3));
        end interpolateInRegion2;
        
      algorithm 
        // Determine upstream density, upstream viscosity, and lambda2
      //  d       := if dp >= 0 then d_a else d_b;
      //  eta     := if dp >= 0 then eta_a else eta_b;
        d       := d_up;
        eta     := eta_up;
        lambda2 := abs(dp)*2*diameter^3*d/(length*eta*eta);
        
        // Determine Re under the assumption of laminar flow
        Re := lambda2/64;
        
        // Modify Re, if turbulent flow
        if Re > Re1 then
           Re :=-2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
           if Re < Re2 then
              Re := interpolateInRegion2(Re, Re1, Re2, Delta, lambda2);
           end if;
        end if;
        
        // Determine mass flow rate
        m_flow := (pi*diameter/4)*eta*(if dp >= 0 then Re else -Re);
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        import Modelica.Math;
                annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi = Modelica.Constants.pi;
        Real Delta = roughness/diameter "Relative roughness";
        SI.ReynoldsNumber Re1 = 745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta) 
          "Re leaving laminar curve";
        SI.ReynoldsNumber Re2 = 4000 "Re entering turbulent curve";
        SI.DynamicViscosity eta "Upstream viscosity";
        SI.Density d "Upstream density";
        SI.ReynoldsNumber Re "Reynolds number";
        Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
        
        function interpolateInRegion2 
           input SI.ReynoldsNumber Re;
           input SI.ReynoldsNumber Re1;
           input SI.ReynoldsNumber Re2;
           input Real Delta;
           output Real lambda2;
           annotation(smoothOrder=1);
          // point lg(lambda2(Re1)) with derivative at lg(Re1)
        protected 
          Real x1 = Math.log10(Re1);
          Real y1 = Math.log10(64*Re1);
          Real yd1=1;
          
          // Point lg(lambda2(Re2)) with derivative at lg(Re2)
          Real aux1=(0.5/Math.log(10))*5.74*0.9;
          Real aux2=Delta/3.7 + 5.74/Re2^0.9;
          Real aux3=Math.log10(aux2);
          Real L2=0.25*(Re2/aux3)^2;
          Real aux4=2.51/sqrt(L2) + 0.27*Delta;
          Real aux5=-2*sqrt(L2)*Math.log10(aux4);
          Real x2 =  Math.log10(Re2);
          Real y2 =  Math.log10(L2);
          Real yd2 = 2 + 4*aux1/(aux2*aux3*(Re2)^0.9);
          
          // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
          Real diff_x=x2 - x1;
          Real m=(y2 - y1)/diff_x;
          Real c2=(3*m - 2*yd1 - yd2)/diff_x;
          Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
          Real dx;
        algorithm 
           dx := Math.log10(Re/Re1);
           lambda2 := 64*Re1*(Re/Re1)^(1 + dx*(c2 + dx*c3));
        end interpolateInRegion2;
      algorithm 
        // Determine upstream density and upstream viscosity
        // d       :=if m_flow >= 0 then d_a else d_b;
        // eta     :=if m_flow >= 0 then eta_a else eta_b;
        d:=d_up;
        eta:=eta_up;
        
        // Determine Re, lambda2 and pressure drop
        Re :=(4/pi)*abs(m_flow)/(diameter*eta);
        lambda2 := if Re <= Re1 then 64*Re else 
                  (if Re >= Re2 then 0.25*(Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else 
                   interpolateInRegion2(Re, Re1, Re2, Delta));
        dp :=length*eta*eta/(2*d*diameter*diameter*diameter)*
             (if m_flow >= 0 then lambda2 else -lambda2);
      end pressureLoss_m_flow;
    end Detailed;
  end WallFriction;
  
  package Systems "Test systems" 
    
    model AlphaSystem_01 
      FluidConcept.BoundaryConditions.AlphaFlowSourceB source(
        T=300,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false) 
                              annotation (extent=[-100,0; -80,20]);
      replaceable model PipeModel = 
          FluidConcept.DistributedFlow.AlphaDistributedFlow;
      PipeModel pipe1(
        diameter=0.01,
        T0=400,
        length=1,
        n=2,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.25e5,
        p0_b=1.22e5) 
        annotation (extent=[-72,0; -52,20]);
      annotation (Diagram);
      FluidConcept.BoundaryConditions.AlphaFlowSourceB sink(
        useMassFlowRate=false,
        p=1e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        usePressureInput=true,
        T=400) 
        annotation (extent=[80,0; 60,20]);
      PipeModel pipe2(
        T0=400,
        length=1,
        diameter=0.005,
        n=2,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.20e5,
        p0_b=1.16e5) 
        annotation (extent=[-20,14; 0,34]);
      PipeModel pipe3(
        diameter=0.01,
        T0=400,
        length=1,
        n=2,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.20e5,
        p0_b=1.16e5) 
        annotation (extent=[-18,-18; 2,2]);
      PipeModel pipe4(
        diameter=0.01,
        T0=400,
        length=1,
        n=2,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.15e5,
        p0_b=1.01e5) 
        annotation (extent=[26,-2; 46,18]);
      Junctions.AlphaJunction alphaJunction(n_b=2,
        T0=400,
        p0=1.21e5,
        V=0.001) 
        annotation (extent=[-46,-2; -26,18]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[66,52; 86,72]);
    equation 
      connect(source.port, pipe1.port_a) annotation (points=[-79.1,9.7; -72.9,
            9.7], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=69,
          rgbfillColor={0,128,255},
          fillPattern=1));
      connect(pipe2.port_b, pipe4.port_a) annotation (points=[0.9,23.7; 16,23.7; 
            16,7.7; 25.1,7.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe3.port_b, pipe4.port_a) annotation (points=[2.9,-8.3; 16,-8.3; 
            16,7.7; 25.1,7.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe4.port_b, sink.port) annotation (points=[46.9,7.7; 54,10;
            59.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe1.port_b, alphaJunction.port_a[1]) annotation (points=[-51.1,
            9.7; -46.55,9.7; -46.55,7.7; -42.1,7.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(alphaJunction.port_b[1], pipe2.port_a) annotation (points=[-30.3,
            6.95; -30.3,14.85; -20.9,14.85; -20.9,23.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(alphaJunction.port_b[2], pipe3.port_a) annotation (points=[-30.3,
            8.45; -30.3,-1.15; -18.9,-1.15; -18.9,-8.3], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(ramp.y, sink.p_in) annotation (points=[87,62; 98,62; 98,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
    end AlphaSystem_01;
    
    model BetaSystem_01 
      
      FluidConcept.BoundaryConditions.BetaFlowSourceB source(
        T=300,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false) 
                              annotation (extent=[-100,0; -80,20]);
      model PipeModel=FluidConcept.DistributedFlow.GammaDistributedFlow;
      FluidConcept.DistributedFlow.BetaDistributedFlow pipe1(
        diameter=0.01,
        T0=400,
        length=1,
        n=2,
        p0_a=1.2e5,
        p0_b=1.16e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed) 
        annotation (extent=[-72,0; -52,20]);
      annotation (Diagram);
      FluidConcept.BoundaryConditions.BetaFlowSourceA sink(
        useMassFlowRate=false,
        p=1e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        usePressureInput=true,
        T=400) 
        annotation (extent=[80,0; 60,20]);
      FluidConcept.DistributedFlow.BetaDistributedFlow pipe2(
        T0=400,
        length=1,
        diameter=0.005,
        n=2,
        p0_b=1.01e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        p0_a=1.15e5,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed) 
        annotation (extent=[-20,14; 0,34]);
      FluidConcept.DistributedFlow.BetaDistributedFlow pipe3(
        diameter=0.01,
        T0=400,
        length=1,
        n=2,
        p0_b=1.01e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        p0_a=1.15e5,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed) 
        annotation (extent=[-18,-18; 2,2]);
      FluidConcept.DistributedFlow.BetaDistributedFlow pipe4(
        diameter=0.01,
        T0=400,
        length=1,
        n=2,
        p0_a=1.2e5,
        p0_b=1.16e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed) 
        annotation (extent=[32,0; 52,20]);
      Junctions.BetaJunction junction1(n_b=2, V=0.00001,
        p0=1.20e5,
        T0=400,
        initType=Modelica_Fluid.Types.Init.InitialValues) 
        annotation (extent=[-48,0; -28,20]);
      Junctions.BetaJunction junction2(n_a=2, n_b=1,
        V=0.00001,
        T0=400,
        initType=Modelica_Fluid.Types.Init.InitialValues) 
        annotation (extent=[4,0; 24,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[56,42; 76,62]);
    equation 
      connect(source.port, pipe1.port_a) annotation (points=[-79.1,9.7; -72.9,
            9.7], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=69,
          rgbfillColor={0,128,255},
          fillPattern=1));
      connect(pipe1.port_b, junction1.port_a[1]) annotation (points=[-51.1,9.7;
            -44.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(junction1.port_b[1], pipe2.port_a) annotation (points=[-32.3,8.95;
            -32.3,16.85; -20.9,16.85; -20.9,23.7], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(junction1.port_b[2], pipe3.port_a) annotation (points=[-32.3,
            10.45; -32.3,0.85; -18.9,0.85; -18.9,-8.3], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe4.port_b, sink.port) annotation (points=[52.9,9.7; 59.1,9.7],
          style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(junction2.port_b[1], pipe4.port_a) annotation (points=[19.7,9.7;
            31.1,9.7],                       style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(junction2.port_a[2], pipe3.port_b) annotation (points=[7.9,10.45; 
            0,10.45; 0,-2; 6,-2; 6,-8.3; 2.9,-8.3], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe2.port_b, junction2.port_a[1]) annotation (points=[0.9,23.7; 
            0.9,16.85; 7.9,16.85; 7.9,8.95], style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(ramp.y, sink.p_in) annotation (points=[77,52; 92,52; 92,50; 96,50;
            96,17.5; 76.8,17.5], style(color=74, rgbcolor={0,0,127}));
    end BetaSystem_01;
    
    model GammaSystem_01 
      
      FluidConcept.BoundaryConditions.GammaFlowSourceB source(
        T=300,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false) 
                              annotation (extent=[-100,0; -80,20]);
      model PipeModel=FluidConcept.DistributedFlow.GammaDistributedFlow;
      FluidConcept.DistributedFlow.GammaDistributedFlow pipe1(
        diameter=0.01,
        T0=400,
        length=1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        p0_a=1.25e5,
        p0_b=1.20e5,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        n=2) 
        annotation (extent=[-72,0; -52,20]);
      annotation (Diagram);
      FluidConcept.BoundaryConditions.GammaFlowSourceA sink(
        useMassFlowRate=false,
        p=1e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        usePressureInput=true,
        T=400) 
        annotation (extent=[80,0; 60,20]);
      FluidConcept.DistributedFlow.GammaDistributedFlow pipe2(
        T0=400,
        length=1,
        diameter=0.005,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        p0_b=1.16e5,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.19e5,
        n=2) 
        annotation (extent=[-20,14; 0,34]);
      FluidConcept.DistributedFlow.GammaDistributedFlow pipe3(
        diameter=0.01,
        T0=400,
        length=1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        p0_b=1.16e5,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.19e5,
        n=2) 
        annotation (extent=[-18,-18; 2,2]);
      FluidConcept.DistributedFlow.GammaDistributedFlow pipe4(
        diameter=0.01,
        T0=400,
        length=1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        p0_a=1.15e5,
        p0_b=1.01e5,
        n=10,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed) 
        annotation (extent=[32,0; 52,20]);
      Junctions.GammaJunction junction(
        n_b=2,
        p0=1.20e5,
        T0=400,
        V=0.001) 
                annotation (extent=[-50,0; -30,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[46,32; 66,52]);
    equation 
      connect(source.port, pipe1.port_a) annotation (points=[-79.1,9.7; -72.9,
            9.7], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=69,
          rgbfillColor={0,128,255},
          fillPattern=1));
      connect(pipe4.port_b, sink.port) annotation (points=[52.9,9.7; 59.1,9.7],
          style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe2.port_b, pipe4.port_a) annotation (points=[0.9,23.7; 18,23.7;
            18,9.7; 31.1,9.7], style(color=3, rgbcolor={0,0,255}));
      connect(pipe4.port_a, pipe3.port_b) annotation (points=[31.1,9.7; 18,9.7;
            18,-8.3; 2.9,-8.3], style(color=3, rgbcolor={0,0,255}));
      connect(pipe1.port_b, junction.port_a[1]) annotation (points=[-51.1,9.7;
            -46.1,9.7],                           style(color=3, rgbcolor={0,0,
              255}));
      connect(junction.port_b[1], pipe2.port_a) annotation (points=[-34.3,8.95;
            -34.3,17.85; -20.9,17.85; -20.9,23.7], style(color=3, rgbcolor={0,0,
              255}));
      connect(junction.port_b[2], pipe3.port_a) annotation (points=[-34.3,10.45;
            -34.3,-8.3; -18.9,-8.3], style(color=3, rgbcolor={0,0,255}));
      connect(ramp.y, sink.p_in) annotation (points=[67,42; 88,42; 88,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
    end GammaSystem_01;
    
    model DeltaSystem_01 
      
      FluidConcept.BoundaryConditions.AlphaFlowSourceB source(
        T=300,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false) 
                              annotation (extent=[-100,0; -80,20]);
      
      FluidConcept.DistributedFlow.DeltaDistributedFlow pipe1(
        diameter=0.01,
        T0=400,
        length=1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.25e5,
        p0_b=1.20e5,
        n=2) 
        annotation (extent=[-72,0; -52,20]);
      annotation (Diagram);
      FluidConcept.BoundaryConditions.AlphaFlowSourceB sink(
        useMassFlowRate=false,
        p=1e5,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        T=400,
        usePressureInput=true) 
        annotation (extent=[80,0; 60,20]);
      FluidConcept.DistributedFlow.DeltaDistributedFlow pipe2(
        T0=400,
        length=1,
        diameter=0.005,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.20e5,
        p0_b=1.16e5,
        n=2) 
        annotation (extent=[-20,14; 0,34]);
      FluidConcept.DistributedFlow.DeltaDistributedFlow pipe3(
        diameter=0.01,
        T0=400,
        length=1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.20e5,
        p0_b=1.16e5,
        n=2) 
        annotation (extent=[-18,-18; 2,2]);
      FluidConcept.DistributedFlow.DeltaDistributedFlow pipe4(
        diameter=0.01,
        T0=400,
        length=1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        redeclare package PressureDrop = FluidConcept.WallFriction.Detailed,
        p0_a=1.15e5,
        p0_b=1.01e5,
        n=2) 
        annotation (extent=[32,0; 52,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[56,42; 76,62]);
    equation 
      connect(source.port, pipe1.port_a) annotation (points=[-79.1,9.7; -72.9,
            9.7], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=69,
          rgbfillColor={0,128,255},
          fillPattern=1));
      connect(pipe4.port_b, sink.port) annotation (points=[52.9,9.7; 59.1,9.7],
          style(
          color=3,
          rgbcolor={0,0,255},
          gradient=3,
          fillColor=71,
          rgbfillColor={85,170,255}));
      connect(pipe1.port_b, pipe2.port_a) annotation (points=[-51.1,9.7; -38,
            9.7; -38,23.7; -20.9,23.7], style(color=3, rgbcolor={0,0,255}));
      connect(pipe2.port_b, pipe4.port_a) annotation (points=[0.9,23.7; 18,23.7;
            18,9.7; 31.1,9.7], style(color=3, rgbcolor={0,0,255}));
      connect(pipe1.port_b, pipe3.port_a) annotation (points=[-51.1,9.7; -38,
            9.7; -38,-8.3; -18.9,-8.3], style(color=3, rgbcolor={0,0,255}));
      connect(pipe4.port_a, pipe3.port_b) annotation (points=[31.1,9.7; 18,9.7;
            18,-8.3; 2.9,-8.3], style(color=3, rgbcolor={0,0,255}));
      connect(ramp.y, sink.p_in) annotation (points=[77,52; 92,52; 92,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
    end DeltaSystem_01;
    
    model AlphaSystem_02 
      FluidConcept.BoundaryConditions.AlphaFlowSourceA source(
        T=300,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                              annotation (extent=[-100,0; -80,20]);
      
      FluidConcept.BoundaryConditions.AlphaFlowSourceA sink(
        useMassFlowRate=false,
        p=1e5,
        usePressureInput=true,
        T=400,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[80,0; 60,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[66,52; 86,72]);
      ControlVolumes.AlphaCV alphaCV(
        length=0.1,
        diameter=0.02,
        p0=1.1e5,
        T0=300,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                     annotation (extent=[-50,0; -30,20]);
      FlowModels.AlphaFM alphaFM(
        diameter=0.02,
        p0_a=1.2e5,
        p0_b=1.1e5,
        m_flow0=0.1,
        length=0.2,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                 annotation (extent=[-70,2; -50,22]);
      FlowModels.AlphaFM alphaFM1(
        length=0.1,
        diameter=0.02,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[-24,22; -4,42]);
      FlowModels.AlphaFM alphaFM2(
        length=0.1,
        diameter=0.02,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[-22,-16; -2,4]);
      ControlVolumes.AlphaCV alphaCV1(
        length=0.1,
        diameter=0.02,
        p0=1.05e5,
        T0=300,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                      annotation (extent=[4,2; 24,22]);
      FlowModels.AlphaFM alphaFM3(
        length=0.1,
        diameter=0.02,
        p0_a=1.05e5,
        p0_b=1.0e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[32,2; 52,22]);
    equation 
      connect(ramp.y, sink.p_in) annotation (points=[87,62; 98,62; 98,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
      connect(source.port, alphaFM.port_a) annotation (points=[-79.1,9.7;
            -74.55,9.7; -74.55,12.7; -66.9,12.7],
                                         style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM.port_b, alphaCV.port_a) annotation (points=[-53.1,12.7;
            -50.55,12.7; -50.55,9.7; -46.9,9.7],
                                         style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM1.port_a, alphaCV.port_b) annotation (points=[-20.9,32.7; -32,
            32.7; -32,9.7; -33.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaCV.port_b, alphaFM2.port_a) annotation (points=[-33.1,9.7; -33.1,
            -5.3; -18.9,-5.3], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM1.port_b, alphaCV1.port_a) annotation (points=[-7.1,32.7; 7.1,
            32.7; 7.1,11.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM2.port_b, alphaCV1.port_a) annotation (points=[-5.1,-5.3; 7.1,
            -5.3; 7.1,11.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaCV1.port_b, alphaFM3.port_a) annotation (points=[20.9,11.7;
            28.45,11.7; 28.45,12.7; 35.1,12.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM3.port_b, sink.port) annotation (points=[48.9,12.7; 53.45,12.7;
            53.45,9.7; 59.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      annotation (Diagram);
    end AlphaSystem_02;
    
    model AlphaPVSystem_02 
      FluidConcept.BoundaryConditions.AlphaFlowSourceA source(
        T=300,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false,
        redeclare package Medium = Modelica.Media.Water.StandardWater) 
                              annotation (extent=[-100,0; -80,20]);
      
      FluidConcept.BoundaryConditions.AlphaFlowSourceA sink(
        useMassFlowRate=false,
        p=1e5,
        usePressureInput=true,
        T=400,
        redeclare package Medium = Modelica.Media.Water.StandardWater) 
        annotation (extent=[80,0; 60,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[66,52; 86,72]);
      ControlVolumes.AlphaCVPortVolume alphaCV(
        diameter=0.01,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T0=300,
        p0=1.1e5,
        length=0.01,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                     annotation (extent=[-52,0; -32,20]);
      FlowModels.AlphaFM alphaFM(
        length=1,
        diameter=0.01,
        p0_a=1.2e5,
        p0_b=1.1e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                 annotation (extent=[-72,-2; -52,18]);
      FlowModels.AlphaFM alphaFM1(
        length=1,
        diameter=0.01,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                  annotation (extent=[-24,22; -4,42]);
      FlowModels.AlphaFM alphaFM2(
        length=1,
        diameter=0.01,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                  annotation (extent=[-22,-16; -2,4]);
      ControlVolumes.AlphaCVPortVolume alphaCV1(
        diameter=0.01,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T0=300,
        p0=1.05e5,
        redeclare package Medium = Modelica.Media.Water.StandardWater,
        length=0.01)                  annotation (extent=[4,2; 24,22]);
      FlowModels.AlphaFM alphaFM3(
        length=1,
        diameter=0.01,
        p0_a=1.05e5,
        p0_b=1.0e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                  annotation (extent=[32,2; 52,22]);
    equation 
      connect(ramp.y, sink.p_in) annotation (points=[87,62; 98,62; 98,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
      connect(source.port, alphaFM.port_a) annotation (points=[-79.1,9.7; -74.55,
            9.7; -74.55,8.7; -68.9,8.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM3.port_b, sink.port) annotation (points=[48.9,12.7; 53.45,12.7;
            53.45,9.7; 59.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      annotation (Diagram);
      connect(alphaFM.port_b, alphaCV.port) annotation (points=[-55.1,8.7; 
            -47.55,8.7; -47.55,9.7; -41.9,9.7],   style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaCV.port, alphaFM1.port_a) annotation (points=[-41.9,9.7; 
            -41.9,31.85; -20.9,31.85; -20.9,32.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaCV.port, alphaFM2.port_a) annotation (points=[-41.9,9.7; 
            -38.95,9.7; -38.95,-5.3; -18.9,-5.3],  style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaCV1.port, alphaFM2.port_b) annotation (points=[14.1,11.7;
            15.05,11.7; 15.05,-5.3; -5.1,-5.3], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaFM1.port_b, alphaCV1.port) annotation (points=[-7.1,32.7;
            14.45,32.7; 14.45,11.7; 14.1,11.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(alphaCV1.port, alphaFM3.port_a) annotation (points=[14.1,11.7;
            25.05,11.7; 25.05,12.7; 35.1,12.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
    end AlphaPVSystem_02;
    
    model GammaSystem_02 
      FluidConcept.BoundaryConditions.GammaFlowSourceA source(
        T=300,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                              annotation (extent=[-100,0; -80,20]);
      
      FluidConcept.BoundaryConditions.GammaFlowSourceA sink(
        useMassFlowRate=false,
        p=1e5,
        usePressureInput=true,
        T=400,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[80,0; 60,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[66,52; 86,72]);
      ControlVolumes.GammaCV CV1(
        diameter=0.01,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T0=300,
        p0=1.1e5,
        length=0.01,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                     annotation (extent=[-50,0; -30,20]);
      FlowModels.GammaFM FM1(
        length=1,
        diameter=0.01,
        p0_a=1.2e5,
        p0_b=1.1e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                 annotation (extent=[-72,-2; -52,18]);
      FlowModels.GammaFM FM2(
        length=1,
        diameter=0.01,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[-24,22; -4,42]);
      FlowModels.GammaFM FM3(
        length=1,
        diameter=0.01,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[-22,-16; -2,4]);
      ControlVolumes.GammaCV CV2(
        diameter=0.01,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T0=300,
        p0=1.05e5,
        length=0.01,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                      annotation (extent=[4,2; 24,22]);
      FlowModels.GammaFM FM4(
        length=1,
        diameter=0.01,
        p0_a=1.05e5,
        p0_b=1.0e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[32,2; 52,22]);
    equation 
      connect(ramp.y, sink.p_in) annotation (points=[87,62; 98,62; 98,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
      connect(source.port, FM1.port_a)     annotation (points=[-79.1,9.7; -74.55,
            9.7; -74.55,8.7; -68.9,8.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(FM4.port_b, sink.port)      annotation (points=[48.9,12.7; 53.45,12.7;
            53.45,9.7; 59.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      annotation (Diagram);
      connect(FM1.port_b, CV1.port) annotation (points=[-55.1,8.7; -47.55,8.7;
            -47.55,9.7; -39.9,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(FM2.port_a, CV1.port) annotation (points=[-20.9,32.7; -20.9,21.35;
            -39.9,21.35; -39.9,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(CV1.port, FM3.port_a) annotation (points=[-39.9,9.7; -29.95,9.7;
            -29.95,-5.3; -18.9,-5.3], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(FM2.port_b, CV2.port) annotation (points=[-7.1,32.7; 4.45,32.7;
            4.45,11.7; 14.1,11.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(CV2.port, FM3.port_b) annotation (points=[14.1,11.7; 5.05,11.7;
            5.05,-5.3; -5.1,-5.3], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(CV2.port, FM4.port_a) annotation (points=[14.1,11.7; 25.05,11.7;
            25.05,12.7; 35.1,12.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
    end GammaSystem_02;
    
    model BetaSystem_02 
      FluidConcept.BoundaryConditions.BetaFlowSourceB source(
        T=300,
        m_flow=0.1,
        p=1.2e5,
        useMassFlowRate=false,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                              annotation (extent=[-100,0; -80,20]);
      
      FluidConcept.BoundaryConditions.BetaFlowSourceA sink(
        useMassFlowRate=false,
        p=1e5,
        usePressureInput=true,
        T=400,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[80,0; 60,20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.5e5,
        duration=1,
        offset=1.0e5,
        startTime=1) annotation (extent=[66,52; 86,72]);
      FlowModels.BetaFM FM1(
        length=1,
        diameter=0.01,
        p0_a=1.2e5,
        p0_b=1.1e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                 annotation (extent=[-72,-2; -52,18]);
      FlowModels.BetaFM FM2(
        length=1,
        diameter=0.01,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[-24,28; -4,48]);
      FlowModels.BetaFM FM3(
        length=1,
        diameter=0.01,
        p0_a=1.1e5,
        p0_b=1.05e5,
        m_flow0=0.05,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[-22,-16; -2,4]);
      FlowModels.BetaFM FM4(
        length=1,
        diameter=0.01,
        p0_a=1.05e5,
        p0_b=1.0e5,
        m_flow0=0.1,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                                  annotation (extent=[32,2; 52,22]);
      ControlVolumes.BetaCV CV1(
        diameter=0.01,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T0=300,
        p0=1.1e5,
        length=0.01,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        na=1,
        nb=2)                        annotation (extent=[-48,0; -28,20]);
      ControlVolumes.BetaCV CV2(
        diameter=0.01,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T0=300,
        p0=1.05e5,
        length=0.01,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        na=2,
        nb=1)                         annotation (extent=[4,4; 24,24]);
    equation 
      connect(ramp.y, sink.p_in) annotation (points=[87,62; 98,62; 98,17.5;
            76.8,17.5], style(color=74, rgbcolor={0,0,127}));
      annotation (Diagram);
      connect(source.port, FM1.port_a) annotation (points=[-79.1,9.7; -72.55,
            9.7; -72.55,8.7; -68.9,8.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(FM4.port_b, sink.port) annotation (points=[48.9,12.7; 54.45,12.7;
            54.45,9.7; 59.1,9.7], style(
          color=3,
          rgbcolor={0,0,255},
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(FM1.port_b, CV1.port_a[1]) annotation (points=[-55.1,8.7; -49.55,
            8.7; -49.55,9.7; -44.9,9.7], style(
          color=3, 
          rgbcolor={0,0,255}, 
          thickness=2));
      connect(CV1.port_b[2], FM3.port_a) annotation (points=[-31.1,10.45; -28,
            10.45; -28,-5.3; -18.9,-5.3], style(
          color=3, 
          rgbcolor={0,0,255}, 
          thickness=2));
      connect(FM3.port_b, CV2.port_a[2]) annotation (points=[-5.1,-5.3; 6,-5.3; 
            6,14.45; 7.1,14.45], style(
          color=3, 
          rgbcolor={0,0,255}, 
          thickness=2));
      connect(CV2.port_b[1], FM4.port_a) annotation (points=[20.9,13.7; 27.45,
            13.7; 27.45,12.7; 35.1,12.7], style(
          color=3, 
          rgbcolor={0,0,255}, 
          thickness=2));
      connect(FM2.port_a, CV1.port_b[1]) annotation (points=[-20.9,38.7; -20.9,
            40; -30,40; -30,8.95; -31.1,8.95], style(
          color=3, 
          rgbcolor={0,0,255}, 
          thickness=2));
      connect(FM2.port_b, CV2.port_a[1]) annotation (points=[-7.1,38.7; -7.1,
            39.35; 7.1,39.35; 7.1,12.95], style(
          color=3, 
          rgbcolor={0,0,255}, 
          thickness=2));
    end BetaSystem_02;
  end Systems;
  
  package BoundaryConditions "Sources and sinks, boundary conditions" 
    partial model PartialFlowSource 
      "Ideal flow source that produces a prescribed mass flow with prescribed temperature and mass fraction" 
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
        annotation (choicesAllMatching = true);
      parameter Boolean useTemperature = true 
        "Use temperature instead of enthalpy";
      parameter Boolean useMassFlowRate = true 
        "Use mass flow rate instead of pressure";
      parameter Boolean useFlowRateInput = false 
        "Get the mass flow rate from the input connector" 
        annotation(Evaluate=true);
      parameter Boolean useTemperatureInput = false 
        "Get the temperature from the input connector" 
        annotation(Evaluate = true, Dialog(enable = useTemperature));
      parameter Boolean useEnthalpyInput = false 
        "Get the specific enthalpy from the input connector" 
        annotation(Evaluate = true, Dialog(enable = not useTemperature));
      parameter Boolean usePressureInput = false 
        "Get pressure from input connector" 
        annotation(Evaluate = true, Dialog(enable = not useMassFlowRate));
      parameter Medium.MassFlowRate m_flow 
        "Fixed mass flow rate going out of the fluid port" 
        annotation (Evaluate = true,
                    Dialog(enable = not useFlowRateInput and useMassFlowRate));
      parameter Medium.AbsolutePressure p "Fixed value of pressure" 
          annotation(Evaluate = true, Dialog(enable = not usePressureInput and not useMassFlowRate));
      parameter Medium.Temperature T "Fixed value of temperature" 
        annotation (Evaluate = true,
                    Dialog(enable = not useTemperatureInput and useTemperature));
      parameter Medium.SpecificEnthalpy h "Fixed value of specific enthalpy" 
        annotation (Evaluate = true, Dialog(enable = not useEnthalpyInput and not useTemperature));
      Medium.BaseProperties medium;
      Modelica.Blocks.Interfaces.RealInput m_flow_in(
        redeclare type SignalType = Medium.MassFlowRate) if useFlowRateInput 
        "Prescribed mass flow rate" 
        annotation (extent=[-128,4; -88,44]);
      Modelica.Blocks.Interfaces.RealInput T_in(
        redeclare type SignalType = Medium.Temperature) if useTemperatureInput 
        "Prescribed fluid temperature" 
        annotation (extent=[-128,-48; -88,-8]);
      Modelica.Blocks.Interfaces.RealInput h_in(
        redeclare type SignalType = Medium.SpecificEnthalpy) if useEnthalpyInput 
        "Prescribed fluid specific enthalpy" 
        annotation (extent=[-90,-87; -50,-47]);
       Modelica.Blocks.Interfaces.RealInput p_in(redeclare type SignalType = 
            Medium.AbsolutePressure) if usePressureInput "Prescribed pressure" 
        annotation (extent=[-88,55; -48,95]);
    protected 
      Modelica.Blocks.Interfaces.RealInput m_flow_in_internal(
        redeclare type SignalType = Medium.MassFlowRate) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput T_in_internal(
        redeclare type SignalType = Medium.Temperature) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput h_in_internal(
        redeclare type SignalType = Medium.SpecificEnthalpy) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput p_in_internal(
        redeclare type SignalType = Medium.AbsolutePressure) 
        "Needed to connect to conditional connector";
      
      annotation (defaultComponentName = "massFlowRate",
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[1,1],
          component=[20, 20],
          scale=0),
        Icon(
          Rectangle(extent=[20, 60; 100, -60], style(
              color=0,
              gradient=2,
              fillColor=8)),
          Rectangle(extent=[38, 40; 100, -40], style(
              color=69,
              gradient=2,
              fillColor=69)),
          Ellipse(extent=[-100, 80; 60, -80], style(fillColor=69, rgbfillColor={0,
                  128,255})),
          Text(extent=[-150,160; 150,110],  string="%name"),
          Text(
            extent=[-205,76; -65,44],
            style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255}),
            string="m_flow"),
          Text(
            extent=[-149,-54; -109,-88],
            style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255}),
            string="T"),
          Text(
            extent=[-118,-81; 2,-122],
            style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255}),
            string="h"),
          Text(
            extent=[-122,118; 18,82],
            style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255}),
            string="p")),
        Window(
          x=0.45,
          y=0.01,
          width=0.44,
          height=0.65),
        Diagram);
      
    equation 
      if not usePressureInput then
        p_in_internal = p;
      else
        connect(p_in, p_in_internal);
      end if;
      if not useFlowRateInput then
        m_flow_in_internal = m_flow;
      else
         connect(m_flow_in, m_flow_in_internal);
      end if;
      if not useTemperatureInput or not useTemperature then
        T_in_internal = T;
      else
        connect(T_in, T_in_internal);
      end if;
      if not useEnthalpyInput or useTemperature then
        h_in_internal = h;
      else
        connect(h_in, h_in_internal);
      end if;
      
      if useTemperature then
        medium.T = T_in_internal;
      else
        medium.h = h_in_internal;
      end if;
    end PartialFlowSource;
    
    model AlphaFlowSourceB 
      extends PartialFlowSource;
      
      parameter Modelica_Fluid.Types.SourceFlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"               annotation(Evaluate=true, Dialog(tab="Advanced"));
      Interfaces.AlphaFlowPortB port(
                       m_flow(max=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.OutOfPort then 
                     0 else  +Modelica.Constants.inf,
                              min=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.InToPort then 
                     0 else  -Modelica.Constants.inf)) 
                                    annotation (extent=[94,-18; 124,12]);
      
    equation 
    if useMassFlowRate then
      port.m_flow = -m_flow_in_internal;
    else
      medium.p = p_in_internal;
    end if;
      port.p = medium.p;
      port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
      annotation (Icon(
          Text(
            extent=[-86,18; 48,-16],
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="alpha")));
    end AlphaFlowSourceB;
    
    model GammaFlowSourceA 
      extends PartialFlowSource;
      
      parameter Modelica_Fluid.Types.SourceFlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"               annotation(Evaluate=true, Dialog(tab="Advanced"));
      
     Interfaces.AlphaFlowPortA port(
                       m_flow(max=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.OutOfPort then 
                     0 else  +Modelica.Constants.inf,
                              min=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.InToPort then 
                     0 else  -Modelica.Constants.inf)) 
                                    annotation (extent=[94,-18; 124,12]);
      
    equation 
    if useMassFlowRate then
      port.m_flow = -m_flow_in_internal;
    else
      medium.p = p_in_internal;
    end if;
      port.p = medium.p;
      port.h = medium.h;
        annotation (Icon(
          Text(
            extent=[-86,18; 48,-16],
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="gamma")));
    end GammaFlowSourceA;
    
    model BetaFlowSourceA 
       extends PartialFlowSource;
      
      parameter Modelica_Fluid.Types.SourceFlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"               annotation(Evaluate=true, Dialog(tab="Advanced"));
      
     Interfaces.BetaFlowPortA port(
                       m_flow(max=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.OutOfPort then 
                     0 else  +Modelica.Constants.inf,
                              min=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.InToPort then 
                     0 else  -Modelica.Constants.inf)) 
                                    annotation (extent=[94,-18; 124,12]);
      parameter Boolean portAtB=true annotation(Evaluate=true);
    equation 
    if useMassFlowRate then
      port.m_flow = -m_flow_in_internal;
    else
      medium.p = p_in_internal;
    end if;
      port.p = medium.p;
      port.h_ab = medium.h;
      
      annotation (Icon(
          Text(
            extent=[-88,20; 46,-14],
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="beta")));
    end BetaFlowSourceA;
    
    model GammaFlowSourceB 
      extends PartialFlowSource;
      
      parameter Modelica_Fluid.Types.SourceFlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"               annotation(Evaluate=true, Dialog(tab="Advanced"));
      
     Interfaces.AlphaFlowPortB port(
                       m_flow(max=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.OutOfPort then 
                     0 else  +Modelica.Constants.inf,
                              min=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.InToPort then 
                     0 else  -Modelica.Constants.inf)) 
                                    annotation (extent=[94,-18; 124,12]);
      
    equation 
    if useMassFlowRate then
      port.m_flow = -m_flow_in_internal;
    else
      medium.p = p_in_internal;
    end if;
      port.p = medium.p;
      port.H_flow = if port.m_flow>=0 then port.m_flow*port.h else port.m_flow*medium.h;
      
      annotation (Icon(
          Text(
            extent=[-86,18; 48,-16],
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="gamma")));
    end GammaFlowSourceB;
    
    model BetaFlowSourceB 
       extends PartialFlowSource;
      
      parameter Modelica_Fluid.Types.SourceFlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"               annotation(Evaluate=true, Dialog(tab="Advanced"));
      
     Interfaces.BetaFlowPortB port(
                       m_flow(max=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.OutOfPort then 
                     0 else  +Modelica.Constants.inf,
                              min=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.InToPort then 
                     0 else  -Modelica.Constants.inf)) 
                                    annotation (extent=[94,-18; 124,12]);
      parameter Boolean portAtB=true annotation(Evaluate=true);
    equation 
    if useMassFlowRate then
      port.m_flow = -m_flow_in_internal;
    else
      medium.p = p_in_internal;
    end if;
      port.p = medium.p;
      port.h_ba = medium.h;
      
      annotation (Icon(
          Text(
            extent=[-88,20; 46,-14],
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="beta")));
    end BetaFlowSourceB;
    
    model AlphaFlowSourceA 
      extends PartialFlowSource;
      
      parameter Modelica_Fluid.Types.SourceFlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"               annotation(Evaluate=true, Dialog(tab="Advanced"));
      Interfaces.AlphaFlowPortA port(
                       m_flow(max=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.OutOfPort then 
                     0 else  +Modelica.Constants.inf,
                              min=if flowDirection==Modelica_Fluid.Types.SourceFlowDirection.InToPort then 
                     0 else  -Modelica.Constants.inf)) 
                                    annotation (extent=[94,-18; 124,12]);
      
    equation 
    if useMassFlowRate then
      port.m_flow = -m_flow_in_internal;
    else
      medium.p = p_in_internal;
    end if;
      port.p = medium.p;
      port.h=medium.h;
      annotation (Icon(
          Text(
            extent=[-86,18; 48,-16],
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="alpha")));
    end AlphaFlowSourceA;
  end BoundaryConditions;
  
  package Junctions "Junction models for branching systems" 
    model AlphaJunction 
      
      annotation (Diagram, Icon(
          Rectangle(extent=[-52,50; 48,-50], style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255})),
          Text(
            extent=[-60,-60; 60,-100],
            style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="%name"),
          Text(
            extent=[-68,12; 60,-12],
            string="alpha",
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}))));
      replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium;
      parameter Boolean steadyState=false annotation(Evaluate=true);
      parameter SI.Volume V=0.001;
      parameter Integer n_a=1 "Number of a ports";
      parameter Integer n_b=1 "Number of b ports";
      parameter Real C=10;
      SI.Mass M;
      SI.InternalEnergy U;
      Medium.BaseProperties medium(preferredMediumStates=true);
      FluidConcept.Interfaces.AlphaFlowPortA[n_a] port_a 
                                       annotation (extent=[-76,-18; -46,12]);
      FluidConcept.Interfaces.AlphaFlowPortB[n_b] port_b 
                                       annotation (extent=[42,-18; 72,12]);
      
    //Initialization
      parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.Init.NoInit 
        "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean useInitialTemperature=true   annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0 annotation(Dialog(tab = "Initialization"));
      
      parameter Medium.Temperature T0 annotation(Dialog(tab = "Initialization",enable=useInitialTemperature));
      parameter Medium.SpecificEnthalpy h0=Medium.specificEnthalpy(Medium.setState_pTX(p0, T0)) 
                                                                              annotation(Dialog(tab = "Initialization",enable=not useInitialTemperature));
      
    initial equation 
      if not Medium.singleState then
        if initType == Modelica_Fluid.Types.Init.SteadyState then
          der(medium.p) = 0;
        elseif initType == Modelica_Fluid.Types.Init.InitialValues then
          medium.p = p0;
        end if;
      end if;
      if initType == Modelica_Fluid.Types.Init.SteadyState then
        medium.h = 0;
      elseif initType == Modelica_Fluid.Types.Init.InitialValues then
        medium.h = h0;
      end if;
    equation 
      for i in 1:n_a loop
        port_a[i].p = medium.p;
        port_a[i].H_flow = semiLinear(
              port_a[i].m_flow,
              port_a[i].h,
              medium.h);
      end for;
      for i in 1:n_b loop
        port_b[i].H_flow = semiLinear(
              port_b[i].m_flow,
              port_b[i].h,
              medium.h);
        port_b[i].p - medium.p = C*port_b[i].m_flow;
      end for;
      M = medium.d*V;
      U = medium.u*M;
      if steadyState then
        0 = sum(port_a.H_flow) + sum(port_b.H_flow);
        0 = sum(port_a.m_flow) + sum(port_b.m_flow);
      else
        der(U) = sum(port_a.H_flow) + sum(port_b.H_flow);
        der(M) = sum(port_a.m_flow) + sum(port_b.m_flow);
      end if;
      
    end AlphaJunction;
    
    model GammaJunction 
      
      annotation (Diagram, Icon(
          Rectangle(extent=[-52,50; 48,-50], style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255})),
          Text(
            extent=[-60,-60; 60,-100],
            style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="%name"),
          Text(
            extent=[-68,12; 60,-12],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="gamma")));
      replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium;
      parameter Boolean steadyState=false annotation(Evaluate=true);
      parameter SI.Volume V=0.001;
      parameter Integer n_a=1 "Number of a ports";
      parameter Integer n_b=1 "Number of b ports";
      parameter Real C=10;
      SI.Mass M;
      SI.InternalEnergy U;
      Medium.BaseProperties medium(preferredMediumStates=true);
      FluidConcept.Interfaces.AlphaFlowPortA[n_a] port_a 
                                       annotation (extent=[-76,-18; -46,12]);
      FluidConcept.Interfaces.AlphaFlowPortB[n_b] port_b 
                                       annotation (extent=[42,-18; 72,12]);
      
    //Initialization
      parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.Init.NoInit 
        "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean useInitialTemperature=true   annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0 annotation(Dialog(tab = "Initialization"));
      
      parameter Medium.Temperature T0 annotation(Dialog(tab = "Initialization",enable=useInitialTemperature));
      parameter Medium.SpecificEnthalpy h0=Medium.specificEnthalpy(
          Medium.setState_pTX(p0, T0))                                        annotation(Dialog(tab = "Initialization",enable=not useInitialTemperature));
      
    initial equation 
      if not Medium.singleState then
        if initType == Modelica_Fluid.Types.Init.SteadyState then
          der(medium.p) = 0;
        elseif initType == Modelica_Fluid.Types.Init.InitialValues then
          medium.p = p0;
        end if;
      end if;
      if initType == Modelica_Fluid.Types.Init.SteadyState then
        medium.h = 0;
      elseif initType == Modelica_Fluid.Types.Init.InitialValues then
        medium.h = h0;
      end if;
    equation 
      for i in 1:n_a loop
        port_a[i].p = medium.p;
        port_a[i].h = medium.h;
        // port_a[i].H_flow = semiLinear(
        //   port_a[i].m_flow, 
        //   port_a[i].h, 
        //   medium.h);
      end for;
      for i in 1:n_b loop
        port_b[i].H_flow = noEvent(if port_b[i].m_flow>=0 then port_b[i].m_flow*port_b[i].h else port_b[i].m_flow*medium.h);
        port_b[i].p - medium.p = C*port_b[i].m_flow;
      end for;
      M = medium.d*V;
      U = medium.u*M;
      if steadyState then
        0 = sum(port_a.H_flow) + sum(port_b.H_flow);
        0 = sum(port_a.m_flow) + sum(port_b.m_flow);
      else
        der(U) = sum(port_a.H_flow) + sum(port_b.H_flow);
        der(M) = sum(port_a.m_flow) + sum(port_b.m_flow);
      end if;
      
    end GammaJunction;
    
    model BetaJunction 
      
      annotation (Diagram, Icon(
          Rectangle(extent=[-52,50; 48,-50], style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255})),
          Text(
            extent=[-60,-60; 60,-100],
            style(
              color=3,
              rgbcolor={0,0,255},
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="%name"),
          Text(
            extent=[-68,12; 60,-12],
            style(
              color=0,
              rgbcolor={0,0,0},
              pattern=0,
              gradient=3,
              fillColor=71,
              rgbfillColor={85,170,255}),
            string="beta")));
      replaceable package Medium = Modelica.Media.Air.DryAirNasa extends 
        Modelica.Media.Interfaces.PartialMedium;
      parameter Boolean steadyState=false annotation(Evaluate=true);
      parameter SI.Volume V=0.001;
      parameter Integer n_a=1 "Number of a ports";
      parameter Integer n_b=1 "Number of b ports";
      parameter Real C=10;
      SI.Mass M;
      SI.InternalEnergy U;
      Medium.BaseProperties medium(preferredMediumStates=true);
      FluidConcept.Interfaces.BetaFlowPortA[n_a] port_a 
                                       annotation (extent=[-76,-18; -46,12]);
      FluidConcept.Interfaces.BetaFlowPortB[n_b] port_b 
                                       annotation (extent=[42,-18; 72,12]);
      
    //Initialization
      parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.Init.NoInit 
        "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean useInitialTemperature=true   annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p0 annotation(Dialog(tab = "Initialization"));
      
      parameter Medium.Temperature T0 annotation(Dialog(tab = "Initialization",enable=useInitialTemperature));
      parameter Medium.SpecificEnthalpy h0=Medium.specificEnthalpy(
          Medium.setState_pTX(p0, T0))                                        annotation(Dialog(tab = "Initialization",enable=not useInitialTemperature));
    protected 
      SI.EnthalpyFlowRate[n_a] H_flow_a;
      SI.EnthalpyFlowRate[n_b] H_flow_b;
      
    initial equation 
      if not Medium.singleState then
        if initType == Modelica_Fluid.Types.Init.SteadyState then
          der(medium.p) = 0;
        elseif initType == Modelica_Fluid.Types.Init.InitialValues then
          medium.p = p0;
        end if;
      end if;
      if initType == Modelica_Fluid.Types.Init.SteadyState then
        medium.h = 0;
      elseif initType == Modelica_Fluid.Types.Init.InitialValues then
        medium.h = h0;
      end if;
    equation 
      for i in 1:n_a loop
        port_a[i].p = medium.p;
        port_a[i].h_ab = medium.h;
        H_flow_a[i] = if noEvent(port_a[i].m_flow>=0) then port_a[i].m_flow*port_a[i].h_ba else port_a[i].m_flow*port_a[i].h_ab;
      end for;
      for i in 1:n_b loop
        port_b[i].h_ba = medium.h;
        port_b[i].p - medium.p = C*port_b[i].m_flow;
        H_flow_b[i] =if  noEvent(port_b[i].m_flow>=0) then port_b[i].m_flow*port_b[i].h_ab else port_b[i].m_flow*port_b[i].h_ba;
      end for;
      
      M = medium.d*V;
      U = medium.u*M;
      if steadyState then
       0 = sum(H_flow_a) + sum(H_flow_b);
        0 = sum(port_a.m_flow) + sum(port_b.m_flow);
      else
        der(U) = sum(H_flow_a) + sum(H_flow_b);
        der(M) = sum(port_a.m_flow) + sum(port_b.m_flow);
      end if;
      
    end BetaJunction;
  end Junctions;
end FluidConcept;
