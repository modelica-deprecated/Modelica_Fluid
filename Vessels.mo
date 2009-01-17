within Modelica_Fluid;
package Vessels "Devices for storing fluid"
   extends Modelica_Fluid.Icons.VariantLibrary;

    model ClosedVolume
    "Volume of fixed size, closed to the ambient, with inlet/outlet ports"
    import Modelica.Constants.pi;

      // Mass and energy balance, ports
      extends Modelica_Fluid.Vessels.BaseClasses.PartialLumpedVessel(
        final fluidVolume = V,
        heatTransfer(surfaceAreas={4*pi*(3/4*V/pi)^(2/3)}));

      parameter SI.Volume V "Volume";

    equation
      Wb_flow = 0;
      for i in 1:nPorts loop
        ports_p_static[i] = medium.p;
      end for;

      annotation (defaultComponentName="volume",
        Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={170,213,255}), Text(
            extent={{-150,12},{150,-18}},
            lineColor={0,0,0},
            textString="V=%V")}),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. 
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. 
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. 
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}),
              graphics));
    end ClosedVolume;

model SimpleTank "Simple tank with inlet/outlet ports"
    import Modelica.Constants.pi;

  // Tank properties
  SI.Height level(stateSelect=StateSelect.prefer, start=level_start)
      "Level height of tank";
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";

  // Mass and energy balance, ports
  extends Modelica_Fluid.Vessels.BaseClasses.PartialLumpedVessel(
    final fluidVolume = V,
    final fluidLevel = level,
    heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
    final initialize_p = false,
    final p_start = p_ambient,
    final use_d_nominal = false,
    final d_nominal = 0);

  // Tank geometry
  parameter SI.Height height "Height of tank";
  parameter SI.Area crossArea "Area of tank";

  // Initialization
  parameter SI.Height level_start(min=0) = 0.5*height
      "Start value of tank level" 
    annotation(Dialog(tab="Initialization"));

  // Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));

equation
  // Total quantities
  V = crossArea*level "Volume of fluid";
  medium.p = p_ambient;

  // Source termsEnergy balance
  if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
    Wb_flow = 0
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
  else
    Wb_flow = -p_ambient*der(V);
  end if;

  assert(level <= height, "Tank is full (level = height = " + String(level) + ")");

  //Determine port properties
  for i in 1:nPorts loop
    ports_p_static[i] = max(0, level - portsData_height[i])*system.g*medium.d + p_ambient;
  end for;

initial equation
  if massDynamics == Types.Dynamics.FixedInitial then
    level = level_start;
  elseif massDynamics == Types.Dynamics.SteadyStateInitial then
    der(level) = 0;
  end if;

    annotation (defaultComponentName="tank",
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-100,100},{100,0}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent=DynamicSelect({{-100,-100},{100,10}}, {{-100,-100},{100,(-100
                 + 200*level/height)}}),
            lineColor={0,127,255},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                0}),
          Text(
            extent={{-129,40},{130,26}},
            lineColor={0,0,0},
            textString="%level_start m"),
          Text(
            extent={{-95,30},{95,5}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level, 
                1, 
                integer(precision)))),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot),
          Text(
            extent={{-126,81},{133,67}},
            lineColor={0,0,0},
            textString="level_start =")}),
      Documentation(info="<HTML>
<p> 
Model of a tank that is open to the ambient at the fixed pressure
<tt>p_ambient</tt>.
</p>
<p>
The vector of connectors <b>ports</b> represents fluid ports at configurable heights, relative to the bottom of tank. 
Fluid can flow either out of or in to each port.
</p>
The following assumptions are made:
<ul>
<li>The tank is filled with a single or multiple-substance medium having a density higher than the density of the ambient medium.</li>
<li>The fluid has uniform density, temperature and mass fractions</li>
<li>No liquid is leaving the tank through the open top; the simulation breaks with an assertion if the liquid level growths over the height.</li>
</ul>
<p>
The port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe. 
The hydraulic resistances <tt>zetas_in</tt> and <tt>zetas_out</tt> determine the dissipative pressure drop between tank and port depending on 
the direction of mass flow. The default values (zetas_in=1, zetas_out=0) assume no dissipation at the tank outlet (ideal smooth opening) and 
total dissipation of kinetic energy at the tank inlet. Larger values are found for sharp edged openings and non-uniform velocity distributions 
in the pipe. A large selection of possible cases are listed in <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. 
</p>
<p>
With the setting <tt>use_portsData=false</tt>, the port pressure represents the static head 
at the height of the respective port. 
The relationship between pressure drop and mass flow rate at the port must then be provided by connected components; 
Heights of ports as well as kinetic and potential energy of fluid enering or leaving are not taken into account anymore. 
</p>   
</HTML>", revisions="<html>
<ul>
<li><i>Dec. 12, 2008</i> by Ruediger Franke: move port definitions 
   to BaseClasses.PartialLumpedVessel; also use energy and mass balance from common base class</li>
<li><i>Dec. 8, 2008</i> by Michael Wetter (LBNL):<br>
Implemented trace substances.</li>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model.</li>
<li><i>Oct. 29, 2007</i> by Carsten Heinrich (ILK Dresden):<br>
Adapted to the new fluid library interfaces: 
<ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
    <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics),
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")));

end SimpleTank;

model TankWithTopPorts
    "Tank with inlet/outlet ports and with inlet ports at the top"

    import Modelica.Constants;
    import Modelica_Fluid.Fittings.BaseClasses.lossConstant_D_zeta;
    import Modelica_Fluid.Utilities.regRoot2;
    import Modelica_Fluid.Vessels.BaseClasses.VesselPortsData;

  SI.Height level(stateSelect=StateSelect.prefer, start=level_start)
      "Fluid level in the tank";

  //Mass and energy balance
  extends Modelica_Fluid.Interfaces.PartialLumpedVolume(
    final fluidVolume = V,
    final initialize_p = false,
    final p_start = p_ambient,
    final use_d_nominal = false,
    final d_nominal = 0);

  //Tank geometry
  parameter SI.Height height "Maximum level of tank before it overflows";
  parameter SI.Area crossArea "Area of tank";
  parameter SI.Volume V0=0 "Volume of the liquid when level = 0";

  //Port definitions
  parameter Integer nTopPorts = 0 "Number of inlet ports above height (>= 1)" 
                                                annotation(Dialog(__Dymola_connectorSizing=true));

  Modelica_Fluid.Interfaces.FluidPorts_a topPorts[nTopPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0, each min=0))
      "Inlet ports over height at top of tank (fluid flows only from the port in to the tank)"
    annotation (Placement(transformation(
        extent={{-20,0},{20,10}},
        origin={0,100})));
/*
    annotation (Placement(transformation(
        extent={{0,-20},{10,20}},
        rotation=90,
        origin={0,100})));
*/

  parameter Integer nPorts = 0
      "Number of inlet/outlet ports (on bottom and on the side)" 
     annotation(Dialog(__Dymola_connectorSizing=true));
  parameter Modelica_Fluid.Vessels.BaseClasses.VesselPortsData portsData[
                                                                      nPorts]
      "Data of inlet/outlet ports at side and bottom of tank";

  Modelica_Fluid.Interfaces.FluidPorts_b ports[nPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0))
      "inlet/outlet ports at bottom or side of tank (fluid flows in to or out of port; a port might be above the fluid level)"
    annotation (Placement(transformation(
        extent={{-20,0},{20,-10}},
        origin={0,-100})));
/*
    annotation (Placement(transformation(
        extent={{0,-20},{-10,20}},
        rotation=90,
        origin={0,-100})));
*/

  //Initialization
  parameter SI.Height level_start(min=0) = 0.5*height
      "Start value of tank level" 
    annotation(Dialog(tab="Initialization"));

  // Heat transfer through boundary
  parameter Boolean use_HeatTransfer = false
      "= true to use the HeatTransfer model" 
      annotation (Dialog(tab="Assumptions", group="Heat transfer"));
  replaceable model HeatTransfer = 
      Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer 
    constrainedby
      Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
      "Wall heat transfer" 
      annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);
  HeatTransfer heatTransfer(
    redeclare final package Medium = Medium,
    final n=1,
    final states = {medium.state},
    surfaceAreas={crossArea+2*sqrt(crossArea*Modelica.Constants.pi)*level},
    final use_k = use_HeatTransfer) 
      annotation (Placement(transformation(
        extent={{-10,-10},{30,30}},
        rotation=90,
        origin={-50,-10})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer 
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  // Advanced
  parameter Real hysteresisFactor(min=0) = 0.1
      "Hysteresis for empty pipe = diameter*hysteresisFactor" 
    annotation(Dialog(tab="Advanced", group="Port properties"));
  parameter Boolean stiffCharacteristicForEmptyPort = false
      "=true, if steep pressure loss characteristic for empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Port properties"), Evaluate=true);
  parameter Real zetaLarge(min=0) = 1e5
      "Large pressure loss factor if mass flows out of empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));
  parameter SI.MassFlowRate m_flow_small(min=0) = 0.01
      "Regularization range at zero mass flow rate" 
    annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));

  //Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));

  // Tank properties
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
  Medium.EnthalpyFlowRate H_flow_top[nTopPorts]
      "Enthalpy flow rates from the top ports in to the tank";
  Medium.EnthalpyFlowRate port_b_H_flow_bottom[nPorts]
      "Enthalpy flow rates from the bottom ports in to the tank";
  Medium.MassFlowRate mXi_flow_top[nTopPorts, Medium.nXi]
      "Substance mass flow rates from the top ports into the tank";
  Medium.MassFlowRate port_b_mXi_flow_bottom[nPorts, Medium.nXi]
      "Substance mass flow rates from the bottom ports into the tank";
  Medium.MassFlowRate mC_flow_top[nTopPorts, Medium.nC]
      "Trace substance mass flow rates from the top ports into the tank";
  Medium.MassFlowRate port_b_mC_flow_bottom[nPorts, Medium.nC]
      "Trace substance mass flow rates from the bottom ports into the tank";
  protected
    SI.Area bottomArea[nPorts];
    SI.Diameter ports_emptyPipeHysteresis[nPorts];
    SI.Length levelAbovePort[nPorts] "Height of fluid over bottom ports";
    Boolean ports_m_flow_out[nPorts](each start = true, each fixed=true);
    Boolean aboveLevel[nPorts] "= true, if level >= ports[i].height";
    Real zetas_out[nPorts];
    Modelica.Blocks.Interfaces.RealInput portsData_diameter[nPorts] = portsData.diameter if nPorts > 0;
    Modelica.Blocks.Interfaces.RealInput portsData_diameter2[nPorts];
    Modelica.Blocks.Interfaces.RealInput portsData_height[nPorts] = portsData.height if nPorts > 0;
    Modelica.Blocks.Interfaces.RealInput portsData_height2[nPorts];
equation
  assert(level <= height, "Tank starts to overflow (level = height = " + String(level) + ")");
  assert(m>=0, "Mass in tank is zero");

  // Compute constant data
  connect(portsData_diameter, portsData_diameter2);
  connect(portsData_height,portsData_height2);

  for i in 1:nPorts loop
      bottomArea[i]=Constants.pi*(portsData_diameter2[i]/2)^2;
      ports_emptyPipeHysteresis[i] = portsData_diameter2[i]*hysteresisFactor;
  end for;

  // Only one connection allowed to a port to avoid unwanted ideal mixing
/*
for i in 1:nTopPorts loop
  assert(cardinality(topPorts[i]) <= 1,"
topPorts[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
 
for i in 1:nPorts loop
  assert(cardinality(ports[i]) <= 1,"
ports[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
*/

  // Total quantities
  medium.p = p_ambient;
  V = crossArea*level + V0 "Volume of fluid";

  // Mass balances
  mb_flow = sum(topPorts.m_flow) + sum(ports.m_flow);
  for i in 1:Medium.nXi loop
    mbXi_flow[i] = sum(mXi_flow_top[:,i]) + sum(port_b_mXi_flow_bottom[:,i]);
  end for;
  for i in 1:Medium.nC loop
    mbC_flow[i]  = sum(mC_flow_top[:,i])  + sum(port_b_mC_flow_bottom[:,i]);
  end for;

  // Energy balance
  Hb_flow = sum(H_flow_top) + sum(port_b_H_flow_bottom);
  Qb_flow = heatTransfer.Q_flows[1];
  if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
    Wb_flow = 0
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
  else
    Wb_flow = -p_ambient*der(V);
  end if;

  // Properties at top ports
    for i in 1:nTopPorts loop
       // It is assumed that fluid flows only from one of the top ports in to the tank and never vice versa
       H_flow_top[i]     = topPorts[i].m_flow*actualStream(topPorts[i].h_outflow);
       mXi_flow_top[i,:] = topPorts[i].m_flow*actualStream(topPorts[i].Xi_outflow);
       mC_flow_top[i,:]  = topPorts[i].m_flow*actualStream(topPorts[i].C_outflow);
       topPorts[i].p     = p_ambient;
       topPorts[i].h_outflow = h_start;
       topPorts[i].Xi_outflow = X_start[1:Medium.nXi];
       topPorts[i].C_outflow  = C_start;
/*
       assert(topPorts[i].m_flow > -1, "Mass flows out of tank via topPorts[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
*/
    end for;

  // Properties at bottom ports
    for i in 1:nPorts loop
       port_b_H_flow_bottom[i]   = ports[i].m_flow*actualStream(ports[i].h_outflow);
       port_b_mXi_flow_bottom[i,:] = ports[i].m_flow*actualStream(ports[i].Xi_outflow);
       port_b_mC_flow_bottom[i,:]  = ports[i].m_flow*actualStream(ports[i].C_outflow);
       aboveLevel[i] = level >= (portsData_height2[i] + ports_emptyPipeHysteresis[i])
                       or pre(aboveLevel[i]) and level >= (portsData_height2[i] - ports_emptyPipeHysteresis[i]);
       levelAbovePort[i] = if aboveLevel[i] then level - portsData_height2[i] else 0;
       ports[i].h_outflow = medium.h;
       ports[i].Xi_outflow = medium.Xi;
       ports[i].C_outflow  = C;

       if stiffCharacteristicForEmptyPort then
          // If port is above fluid level, use large zeta if fluid flows out of port (= small mass flow rate)
          zetas_out[i] = 1 + (if aboveLevel[i] then 0 else zetaLarge);
          ports[i].p = p_ambient + levelAbovePort[i]*system.g*medium.d
                               + Modelica_Fluid.Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                     lossConstant_D_zeta(portsData_diameter2[i], 0.01)/medium.d,
                                     lossConstant_D_zeta(portsData_diameter2[i], zetas_out[i])/medium.d);
          ports_m_flow_out[i] = false;

       else
          // Handling according to Remelhe/Poschlad
          ports_m_flow_out[i] = (pre(ports_m_flow_out[i]) and not ports[i].p>p_ambient)
                                     or ports[i].m_flow < -1e-6;
         if aboveLevel[i] then
             ports[i].p = p_ambient + levelAbovePort[i]*system.g*medium.d -
                               smooth(2,noEvent(if ports[i].m_flow < 0 then ports[i].m_flow^2/
                                     (2*medium.d*bottomArea[i]^2) else 0));
         else
            if pre(ports_m_flow_out[i]) then
               ports[i].m_flow = 0;
            else
               ports[i].p = p_ambient;
            end if;
         end if;
          zetas_out[i] =0;
       end if;
     end for;

initial equation
    for i in 1:nPorts loop
       pre(aboveLevel[i]) = level_start >= portsData_height2[i];
    end for;

    if massDynamics == Types.Dynamics.FixedInitial then
      level = level_start;
    elseif massDynamics == Types.Dynamics.SteadyStateInitial then
      der(level) = 0;
    end if;

    annotation (defaultComponentName="tank",
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent=DynamicSelect({{-100,-100},{100,0}}, {{-100,-100},{100,(-100
                 + 200*level/height)}}),
            lineColor={0,127,255},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-94,19},{96,-1}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level, 
                1, 
                3))),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot),
          Text(
            extent={{-94,90},{95,60}},
            lineColor={0,0,255},
            textString="%name"),
          Text(
            extent={{-95,-85},{95,-65}},
            lineColor={0,0,0},
            textString="%level_start"),
          Text(
            extent={{-95,-55},{95,-35}},
            lineColor={0,0,0},
            textString="level_start ="),
          Text(
            extent={{-95,50},{95,30}},
            lineColor={0,0,0},
            textString="level ="),
          Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                0})}),
      Documentation(info="<HTML>
<p> 
Model of a tank that is open to the environment at the fixed pressure
<tt>p_ambient</tt>. 
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.
</p> 
 
<p>
At the top of the tank over the maximal fill level <b>height</b> 
a vector of FluidPorts, called <b>topPorts</b>, is present.
The assumption is made that fluid flows always in to the tank via these
ports (and never back in to the connector).
</p>
 
<p>
The vector of connectors <b>ports</b> are fluid ports at the bottom
and side of the tank at a defineable height. Fluid can flow either out
of or in to this port. The fluid level of the tank may be below
one of these ports. This case is approximated by introducing a
large pressure flow coefficient so that the mass flow rate
through this port is very small in this case.
</p>
 
<p>
If the tank starts to over flow (i.e., level > height), an
assertion is triggered.
</p>
 
<p>
When the diagram layer is open in the plot environment, the
level of the tank is dynamically visualized. Note, the speed
of the diagram animation in Dymola can be set via command
<b>animationSpeed</b>(), e.g., animationSpeed(speed = 10)
</p>
</HTML>", revisions="<html>
<ul>
<li><i>Dec. 12, 2008</i> by Ruediger Franke: replace energy and mass balances with
   common definition in BaseClasses.PartialLumpedVolume</li>
<li><i>Dec. 8, 2008</i> by Michael Wetter (LBNL):<br>
Implemented trace substances and missing equation for outflow of multi substance media at top port.</li>
<li><i>Jul. 29, 2006</i> by Martin Otter (DLR):<br> 
   Improved handling of ports that are above the fluid level and
   simpler implementation.</li>
 
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model but with several improvements
   (top, bottom, side ports; correctly treating kinetic energy for outlet
   and total dissipation for inlet; ports can be above the fluid level).</li>
</ul>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics),
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")));
equation

    connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
        points={{-100,0},{-87,0},{-87,8.88178e-016},{-74,8.88178e-016}},
        color={191,0,0},
        smooth=Smooth.None));
end TankWithTopPorts;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

      partial model PartialLumpedVessel
      "Lumped volume with a vector of fluid ports and replaceable heat transfer model"
        extends Modelica_Fluid.Interfaces.PartialLumpedVolume;

        // Port definitions
        parameter Integer nPorts=0 "Number of ports" 
          annotation(Evaluate=true, Dialog(__Dymola_connectorSizing=true, tab="General",group="Ports"));
        Interfaces.FluidPorts_b ports[nPorts](redeclare each package Medium = Medium)
        "Fluid outlets" 
          annotation (Placement(transformation(extent={{-40,-10},{40,10}},
            origin={0,-100})));
      /*
    annotation (Placement(transformation(extent={{-10,-40},{10,40}},
      rotation=-90,
      origin={0,-100})));
*/

        input SI.Height fluidLevel = Modelica.Constants.inf
        "level of fluid in the vessel for treating heights of ports";
        Medium.AbsolutePressure[nPorts] ports_p_static
        "static pressures at the ports, inside the vessel";

        // Port properties
        parameter Boolean use_portsData=true
        "= false to neglect pressure loss and kinetic energy" 
          annotation(Evaluate=true, Dialog(tab="General",group="Ports"));
        parameter Modelica_Fluid.Vessels.BaseClasses.VesselPortsData portsData[nPorts] if use_portsData
        "Data of inlet/outlet ports" 
          annotation(Dialog(tab="General",group="Ports",enable= use_portsData));
        parameter Real[nPorts] zetas_in(min=0, max=1)=fill(0.9, nPorts)
        "Hydraulic resistance into volume, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
          annotation(Dialog(tab="General",group="Ports",enable= use_portsData));
        parameter Real[nPorts] zetas_out(min=0, max=1)=fill(0.5, nPorts)
        "Hydraulic resistance out of volume, 0 for ideal smooth outlet" 
          annotation(Dialog(tab="General",group="Ports",enable= use_portsData));

        parameter SI.MassFlowRate m_flow_small(min=0) = 0.01
        "Regularization range at zero mass flow rate" 
          annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));
      /*
  parameter Medium.AbsolutePressure dp_small = 1 
    "Turbulent flow if |dp| >= dp_small (regularization of zero flow)" 
    annotation(Dialog(tab="Advanced",group="Ports"));
*/
        Medium.EnthalpyFlowRate ports_H_flow[nPorts];
        Medium.MassFlowRate ports_mXi_flow[nPorts,Medium.nXi];
        Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
        "Substance mass flows through ports";
        Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts,Medium.nC];
        Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
        "Trace substance mass flows through ports";

        // Heat transfer through boundary
        parameter Boolean use_HeatTransfer = false
        "= true to use the HeatTransfer model" 
            annotation (Dialog(tab="Assumptions", group="Heat transfer"));
        replaceable model HeatTransfer = 
            Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer 
          constrainedby
        Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
        "Wall heat transfer" 
            annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);
        HeatTransfer heatTransfer(
          redeclare final package Medium = Medium,
          final n=1,
          final states = {medium.state},
          final use_k = use_HeatTransfer) 
            annotation (Placement(transformation(
              extent={{-10,-10},{30,30}},
              rotation=90,
              origin={-50,-10})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer 
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

        // Conservation of kinetic energy
        Medium.Density[nPorts] portDensities
        "densites of the fluid at the device boudary";
        SI.Velocity[nPorts] portVelocities
        "velocities of fluid flow at device boundary";
        SI.EnergyFlowRate[nPorts] ports_E_flow
        "flow of kinetic and potential energy at device boundary";

        Real[nPorts] r
        "curve parameters for port flows vs. port pressure drops";
        Real[nPorts] s
        "curve parameters for port flows vs. stored fluid volume";

        // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
        // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
        // providing portsData_diameter and portsData_height, independend of the use_portsData setting.
        // Note: this moreover serves as work-around if a tool does not support a zero sized portsData record.
    protected
        Modelica.Blocks.Interfaces.RealInput[nPorts]
        portsData_diameter_internal =                                              portsData.diameter if use_portsData and nPorts > 0;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal = portsData.height if use_portsData and nPorts > 0;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;

        SI.Area[nPorts] portAreas = {Modelica.Constants.pi/4*portsData_diameter[i]^2 for i in 1:nPorts};

      equation
        mb_flow = sum(ports.m_flow);
        mbXi_flow = sum_ports_mXi_flow;
        mbC_flow  = sum_ports_mC_flow;
        Hb_flow = sum(ports_H_flow) + sum(ports_E_flow);
        Qb_flow = heatTransfer.Q_flows[1];

        // Only one connection allowed to a port to avoid unwanted ideal mixing
        for i in 1:nPorts loop
          assert(cardinality(ports[i]) <= 1,"
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");
        end for;
        // Check for correct solution
        assert(fluidLevel > -1e-6, "Fluid level is below zero meaning that the solution failed.");

        // Boundary conditions

        // treatment of conditional portsData
        connect(portsData_diameter, portsData_diameter_internal);
        connect(portsData_height, portsData_height_internal);
        if not use_portsData then
          portsData_diameter = zeros(nPorts);
          portsData_height = zeros(nPorts);
        end if;

        // actual definition of port variables
        for i in 1:nPorts loop
          if use_portsData then
            // dp = 0.5*zeta*d*v*|v|
            // Note: assume ports_p_static for portDensities to avoid algebraic loops for ports.p
            portDensities[i] = noEvent(Medium.density(Medium.setState_phX(ports_p_static[i], actualStream(ports[i].h_outflow), actualStream(ports[i].Xi_outflow))));
            portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/portDensities[i]);
          else
            // an infinite port diameter is assumed
            portDensities[i] = medium.d;
            portVelocities[i] = 0;
          end if;
        end for;
        for i in 1:nPorts loop
          if fluidLevel >= portsData_height[i] then
            // regular operation: fluidLevel is above ports[i]
            if use_portsData then
              ports[i].p = ports_p_static[i] + (0.5/portAreas[i]^2*Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                           (1 - zetas_in[i])/portDensities[i],
                                           (1 + zetas_out[i])/medium.d));
              /*
        // alternative formulation m_flow=f(dp); not allowing the ideal zetas_in[i]=1 though
        ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - ports_p_static[i], dp_small,
                                     2*portDensities[i]/(1 - zetas_in[i]),
                                     2*medium.d/(1 + zetas_out[i])));
        */
            else
              ports[i].p = ports_p_static[i];
            end if;
            r[i] = fluidLevel - portsData_height[i];
            s[i] = fluidLevel - portsData_height[i];
          elseif s[i] > 0 or r[i] > 0 then
            // ports[i] is above fluidLevel and has inflow
            ports[i].p = ports_p_static[i];
            r[i] = ports[i].m_flow;
            s[i] = ports[i].m_flow;
          else
            // ports[i] is above fluidLevel, preventing outflow
            ports[i].m_flow = 0;
            r[i] = (ports[i].p - ports_p_static[i])/Medium.p_default*(portsData_height[i] - fluidLevel);
            s[i] = fluidLevel - portsData_height[i];
          end if;

          ports[i].h_outflow  = medium.h;
          ports[i].Xi_outflow = medium.Xi;
          ports[i].C_outflow  = C;

          ports_H_flow[i] = ports[i].m_flow * actualStream(ports[i].h_outflow)
          "Enthalpy flow";
          ports_E_flow[i] = ports[i].m_flow*(0.5*portVelocities[i]*portVelocities[i] + system.g*portsData_height[i])
          "Flow of kinetic and potential energy";
          ports_mXi_flow[i,:] = ports[i].m_flow * actualStream(ports[i].Xi_outflow)
          "Component mass flow";
          ports_mC_flow[i,:]  = ports[i].m_flow * actualStream(ports[i].C_outflow)
          "Trace substance mass flow";
        end for;

        for i in 1:Medium.nXi loop
          sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:,i]);
        end for;

        for i in 1:Medium.nC loop
          sum_ports_mC_flow[i]  = sum(ports_mC_flow[:,i]);
        end for;

        connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
            points={{-100,0},{-87,0},{-87,8.88178e-016},{-74,8.88178e-016}},
            color={191,0,0},
            smooth=Smooth.None));
       annotation (
        Documentation(info="<html>
<p>
This base class extends PartialLumpedVolume with a vector of fluid ports and a replaceable wall HeatTransfer model.
It assumes a perfectly mixed volume without kinetic energy in the fluid, i.e. kinetic energy dissipates into the internal energy.
</p>
Each port has a (hydraulic) diameter and a height above the bottom of the vessel, which can be configured using the <b><tt>portsData</tt></b> record.
Alternatively the impact of port geometries can be neglected with <tt>use_portsData=false</tt>. This might be useful for early
design studies. Note that this means to assume an infinite port diameter at the bottom of the vessel. 
Pressure drops and heights of the ports as well as kinetic and potential energy fluid enering or leaving the vessel are neglected then.
<p>
An extending model should use the predefined variables <b><tt>portsData_diameter[nPorts]</tt></b> and <b><tt>portsData_height[nPorts]</tt></b>,
instead of accessing the <tt>portsData</tt> record, as an access to <tt>portsData</tt> may fail for <tt>use_portsData=false</tt> or <tt>nPorts=0</tt>.
The following variables need to be defined by an extending model:
<ul>
<li><tt>input fluidVolume</tt>, the volume of the fluid in the vessel,</li>
<li><tt>input fluidLevel</tt>, the level the fluid in the vessel, which is needed for the treatment of <tt>portsData_height[nPorts]</tt>, and</li>
<li><tt>Wb_flow</tt>, work term of the energy balance, e.g. p*der(V) if the volume is not constant or stirrer power.</li>
</ul>
</html>",       revisions="<html>
<ul>
<li><i>Jan. 2009</i> by R&uuml;diger Franke: extended with
   <ul><li>portsData record and threat configurable port heights,</li>
       <li>consideration of kinetic energy of fluid entering or leaving in energy balance</li>
   </ul>
</li>
<li><i>Dec. 2008</i> by R&uuml;diger Franke: derived from OpenTank, in order to make general use of configurable port diameters</i>
</ul>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}),
                graphics),
          Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}), graphics={Text(
              extent={{-150,110},{150,150}},
              textString="%name",
              lineColor={0,0,255})}));

      end PartialLumpedVessel;

  package HeatTransfer "HeatTransfer models for vessels"

    partial model PartialVesselHeatTransfer
        "Base class for vessel heat transfer models"
      extends Modelica_Fluid.Interfaces.PartialHeatTransfer;

      annotation(Documentation(info="<html>
Base class for vessel heat transfer models.
</html>"),Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                  {100,100}}), graphics={Ellipse(
                extent={{-60,64},{60,-56}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={232,0,0}), Text(
                extent={{-38,26},{40,-14}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={232,0,0},
                textString="%name")}));
    end PartialVesselHeatTransfer;

    model IdealHeatTransfer
        "IdealHeatTransfer: Ideal heat transfer without thermal resistance"
      extends PartialVesselHeatTransfer;

    equation
      Ts = heatPorts.T;

      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end IdealHeatTransfer;

    model ConstantHeatTransfer
        "ConstantHeatTransfer: Constant heat transfer coefficient"
      extends PartialVesselHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0
          "constant heat transfer coefficient";

    equation
      Q_flows = {(alpha0+k)*surfaceAreas[i]*(heatPorts[i].T - Ts[i]) for i in 1:n};

      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient.
</html>"));
    end ConstantHeatTransfer;
    annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));

  end HeatTransfer;

    record VesselPortsData "Data to describe inlet/outlet ports at vessels"
          extends Modelica.Icons.Record;

      parameter SI.Diameter diameter
        "Inner (hydraulic) diameter of inlet/outlet port";
      parameter SI.Height height = 0 "Height over the bottom of the vessel";
    end VesselPortsData;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Vessels;
