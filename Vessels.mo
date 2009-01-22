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
        vessel_ps_static[i] = medium.p;
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
<p>
Ideally mixed volume of constant size with two fluid ports and one medium model. 
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model if <code>use_portsData=false</code>. 
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. 
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</p>
<p>
If <code>use_portsData=true</code>, the port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe. 
The hydraulic resistances <tt>portsData.zeta_in</tt> and <tt>portsData.zeta_out</tt> determine the dissipative pressure drop between volume and port depending on 
the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume an ideal smooth outlet and dissipation for inlet flow. 
Different values are found for sharp edged openings and non-uniform velocity distributions 
in the pipe. Further information can be found in <a href=\"Modelica://Modelica_Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. 
</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}),
              graphics));
    end ClosedVolume;

model SimpleTank "Simple tank with inlet/outlet ports"
    import Modelica.Constants.pi;

  // Tank properties
  SI.Height level(stateSelect=StateSelect.prefer, start=max(level_start, Modelica.Constants.eps))
      "Level height of tank";
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";

  // Tank geometry
  parameter SI.Height height "Height of tank";
  parameter SI.Area crossArea "Area of tank";

  // Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Assumptions", group = "Ambient"));
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Assumptions", group = "Ambient"));

  // Initialization
  parameter SI.Height level_start(min=0) = 0.5*height
      "Start value of tank level" 
    annotation(Dialog(tab="Initialization"));

  // Mass and energy balance, ports
  extends Modelica_Fluid.Vessels.BaseClasses.PartialLumpedVessel(
    final fluidVolume = V,
    final fluidLevel = level,
    final fluidLevel_max = height,
    heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
    final initialize_p = false,
    final p_start = p_ambient);

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

  //Determine port properties
  for i in 1:nPorts loop
    vessel_ps_static[i] = max(0, level - portsData_height[i])*system.g*medium.d + p_ambient;
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
            extent={{-95,75},{95,55}},
            lineColor={0,0,0},
            textString="level ="),
          Text(
            extent={{-95,40},{95,20}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level, 
                1, 
                3))),
          Text(
            extent={{-95,-40},{95,-20}},
            lineColor={0,0,0},
            textString="level_start ="),
          Text(
            extent={{-95,-75},{95,-55}},
            lineColor={0,0,0},
            textString="%level_start"),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot)}),
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
The hydraulic resistances <tt>portsData.zeta_in</tt> and <tt>portsData.zeta_out</tt> determine the dissipative pressure drop between tank and port depending on 
the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume an ideal smooth outlet and dissipation for inlet flow. 
Different values are found for sharp edged openings and non-uniform velocity distributions 
in the pipe. Further information can be found in <a href=\"Modelica://Modelica_Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. 
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

  //Tank geometry
  parameter SI.Height height "Maximum level of tank before it overflows";
  parameter SI.Area crossArea "Area of tank";
  parameter SI.Volume V0=0 "Volume of the liquid when level = 0";

  //Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Assumptions", group = "Ambient"));
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Assumptions", group = "Ambient"));

  //Initialization
  parameter SI.Height level_start(min=0) = 0.5*height
      "Start value of tank level" 
    annotation(Dialog(tab="Initialization"));

  //Mass and energy balance
  extends Modelica_Fluid.Interfaces.PartialLumpedVolume(
    final fluidVolume = V,
    final initialize_p = false,
    final p_start = p_ambient);

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
  parameter SI.MassFlowRate m_flow_small(min=0) = system.m_flow_small
      "Regularization range at zero mass flow rate" 
    annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));

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

        input SI.Height fluidLevel = 0
        "level of fluid in the vessel for treating heights of ports";
        parameter SI.Height fluidLevel_max = 1
        "maximum level of fluid in the vessel";
        Medium.AbsolutePressure[nPorts] vessel_ps_static
        "static pressures inside the vessel at the height of the corresponding ports, zero flow velocity";

        // Port properties
        parameter Boolean use_portsData=true
        "= false to neglect pressure loss and kinetic energy" 
          annotation(Evaluate=true, Dialog(tab="General",group="Ports"));
        parameter Modelica_Fluid.Vessels.BaseClasses.VesselPortsData[nPorts]
        portsData if   use_portsData "Data of inlet/outlet ports" 
          annotation(Dialog(tab="General",group="Ports",enable= use_portsData));

        parameter SI.MassFlowRate m_flow_small(min=0) = system.m_flow_small
        "Regularization range at zero mass flow rate" 
          annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));
      /*
  parameter Medium.AbsolutePressure dp_small = system.dp_small 
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

        // Note: should use fluidLevel_start - portsData.height
        Real[nPorts] s(each start = fluidLevel_max)
        "curve parameters for port flows vs. port pressures";
        Real[nPorts] ports_penetration
        "penetration of port with fluid, depending on fluid level and port diameter";

        // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
        // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
        // providing portsData_diameter and portsData_height, independend of the use_portsData setting.
        // Note: this moreover serves as work-around if a tool does not support a zero sized portsData record.
    protected
        Modelica.Blocks.Interfaces.RealInput[nPorts]
        portsData_diameter_internal =                                              portsData.diameter if use_portsData and nPorts > 0;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal = portsData.height if use_portsData and nPorts > 0;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in_internal = portsData.zeta_in if use_portsData and nPorts > 0;
        Modelica.Blocks.Interfaces.RealInput[nPorts]
        portsData_zeta_out_internal =                                              portsData.zeta_out if use_portsData and nPorts > 0;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in;
        Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out;

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
        assert(fluidLevel <= fluidLevel_max, "Vessel is overflowing (fluidLevel > fluidLevel_max = " + String(fluidLevel) + ")");
        assert(fluidLevel > -1e-6*fluidLevel_max, "Fluid level (= " + String(fluidLevel) + ") is below zero meaning that the solution failed.");

        // Boundary conditions

        // treatment of conditional portsData
        connect(portsData_diameter, portsData_diameter_internal);
        connect(portsData_height, portsData_height_internal);
        connect(portsData_zeta_in, portsData_zeta_in_internal);
        connect(portsData_zeta_out, portsData_zeta_out_internal);
        if not use_portsData then
          portsData_diameter = zeros(nPorts);
          portsData_height = zeros(nPorts);
          portsData_zeta_in = 1*zeros(nPorts);
          portsData_zeta_out = -1*ones(nPorts);
        end if;

        // actual definition of port variables
        for i in 1:nPorts loop
          if use_portsData then
            // dp = 0.5*zeta*d*v*|v|
            // Note: assume vessel_ps_static for portDensities to avoid algebraic loops for ports.p
            portDensities[i] = noEvent(Medium.density(Medium.setState_phX(vessel_ps_static[i], actualStream(ports[i].h_outflow), actualStream(ports[i].Xi_outflow))));
            portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/portDensities[i]);
            // Note: the penetration should not go too close to zero as this would prevent a vessel from running empty
            ports_penetration[i] = Utilities.regStep(fluidLevel - portsData_height[i] - 0.1*portsData_diameter[i], 1, 1e-3, 0.1*portsData_diameter[i]);
          else
            // an infinite port diameter is assumed
            portDensities[i] = medium.d;
            portVelocities[i] = 0;
            ports_penetration[i] = 1;
          end if;
          // fluid flow through ports
          if fluidLevel >= portsData_height[i] then
            // regular operation: fluidLevel is above ports[i]
            // Note: >= covers default values of zero as well
            if use_portsData then
              /* Without regularization
        ports[i].p = vessel_ps_static[i] + 0.5*ports[i].m_flow^2/portAreas[i]^2 
                      * noEvent(if ports[i].m_flow>0 then (zeta_in[i] - 1)/portDensities[i] else -(1+zeta_out[i])/medium.d);
        */

              ports[i].p = vessel_ps_static[i] + (0.5/portAreas[i]^2*Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                           (portsData_zeta_in[i] - 1)/portDensities[i]*ports_penetration[i],
                                           (1 + portsData_zeta_out[i])/medium.d/ports_penetration[i]));
              /*
        // alternative formulation m_flow=f(dp); not allowing the ideal portsData_zeta_in[i]=1 though
        ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - vessel_ps_static[i], dp_small,
                                     2*portDensities[i]/(portsData_zeta_in[i] - 1),
                                     2*medium.d/(1 + portsData_zeta_out[i])));
        */
            else
              ports[i].p = vessel_ps_static[i];
            end if;
            s[i] = fluidLevel - portsData_height[i];
          elseif s[i] > 0 or portsData_height[i] >= fluidLevel_max then
            // ports[i] is above fluidLevel and has inflow
            ports[i].p = vessel_ps_static[i];
            s[i] = ports[i].m_flow;
          else
            // ports[i] is above fluidLevel, preventing outflow
            ports[i].m_flow = 0;
            s[i] = (ports[i].p - vessel_ps_static[i])/Medium.p_default*(portsData_height[i] - fluidLevel);
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
The following variables need to be defined by an extending model:
<ul>
<li><tt>input fluidVolume</tt>, the volume of the fluid in the vessel,</li>
<li><tt>vessel_ps_static[nPorts]</tt>, the static pressures inside the vessel at the height of the corresponding ports, at zero flow velocity, and</li>
<li><tt>Wb_flow</tt>, work term of the energy balance, e.g. p*der(V) if the volume is not constant or stirrer power.</li>
</ul>
Optionally the fluid level may vary in the vessel, which effects the flow through the ports at configurable <tt>portsData_height[nPorts]</tt>. 
This is why an extending model with varying fluid level needs to define:
<ul>
<li><tt>input fluidLevel</tt>, the level the fluid in the vessel, and</li>
<li><tt>input fluidLevel_max</tt>, the maximum level that must not be exceeded. Ports at or above fluidLevel_max can only receive inflow.</li>
</ul>
<p>
An extending model should not access the <tt>portsData</tt> record defined in the configuration dialog,
as an access to <tt>portsData</tt> may fail for <tt>use_portsData=false</tt> or <tt>nPorts=0</tt>.
Instead the predefined variables 
<ul>
<li><tt>portsData_diameter[nPorts]</tt></li>,
<li><tt>portsData_height[nPorts]</tt></li>,
<li><tt>portsData_zeta_in[nPorts]</tt></li>, and
<li><tt>portsData_zeta_out[nPorts]</tt></li>
</ul>
should be used, if needed by an extending model.
</p>
</html>",       revisions="<html>
<ul>
<li><i>Jan. 2009</i> by R&uuml;diger Franke: extended with
   <ul><li>portsData record and threat configurable port heights,</li>
       <li>consideration of kinetic and potential energy of fluid entering or leaving in energy balance</li>
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

    record VesselPortsData "Data to describe inlet/outlet ports at vessels:
    diameter -- Inner (hydraulic) diameter of inlet/outlet port
    height -- Height over the bottom of the vessel
    zeta_out -- Hydraulic resistance out of vessel, default 0.5 for mounted flush with the wall
    zeta_in -- Hydraulic resistance into vessel, default 1.04 for small port diameter"
          extends Modelica.Icons.Record;
      parameter SI.Diameter diameter
        "Inner (hydraulic) diameter of inlet/outlet port";
      parameter SI.Height height = 0 "Height over the bottom of the vessel";
      parameter Real zeta_out(min=0)=0.5
        "Hydraulic resistance out of vessel, default 0.5 for mounted flush with the wall";
      parameter Real zeta_in(min=0)=1.04
        "Hydraulic resistance into vessel, default 1.04 for small port diameter";
      annotation (preferredView="info", Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Vessel Port Data</font></h3>
<p>
This record describes the <b>ports</b> of a <b>vessel</b>. The variables in it are mostly self-explanatory (see list below); only the &zeta; loss factors <code>zeta_inlet</code> and <code>zeta_outlet</code> are discussed further. All data is quoted from Idelchik (1994).
</p>
 
<h4><font color=\"#008000\">Outlet Coefficients</font></h4>
 
<p>
If a <b>straight pipe with constant cross section is mounted flush with the wall</b>, its outlet pressure loss coefficient will be <code>zeta_out[i] = 0.5</code> (Idelchik, p. 160, Diagram 3-1, paragraph 2).
</p>
<p>
If a <b>straight pipe with constant cross section is mounted into a vessel such that the entrance into it is at a distance</b> <code>b</code> from the wall (inside) the following table can be used. Herein, &delta; is the tube wall thickness (Idelchik, p. 160, Diagram 3-1, paragraph 1).
</p> 
 
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for outlets, entrance at a distance from wall</caption>
  <tr>
    <td></td> <td>   </td><th colspan=\"5\" align=\"center\"> b / D_hyd  </th>
  </tr>
  <tr>
    <td></td> <td>   </td><th> 0.000 </th><th> 0.005 </th><th> 0.020 </th><th> 0.100 </th><th> 0.500-&#8734; </th>
  </tr>
  <tr>
     <th rowspan=\"5\" valign=\"middle\">&delta; / D_hyd</th> <th> 0.000 </th><td> 0.50 </td><td> 0.63  </td><td> 0.73  </td><td> 0.86  </td><td>      1.00     </td>
  </tr>
  <tr>
              <th> 0.008 </th><td> 0.50 </td><td> 0.55  </td><td> 0.62  </td><td> 0.74  </td><td>      0.88     </td>
  </tr>
  <tr>
              <th> 0.016 </th><td> 0.50 </td><td> 0.51  </td><td> 0.55  </td><td> 0.64  </td><td>      0.77     </td>
  </tr>
  <tr>
              <th> 0.024 </th><td> 0.50 </td><td> 0.50  </td><td> 0.52  </td><td> 0.58  </td><td>      0.68     </td>
  </tr>
  <tr>
              <th> 0.040 </th><td> 0.50 </td><td> 0.50  </td><td> 0.51  </td><td> 0.51  </td><td>      0.54     </td>
  </tr>
</table>
 
<p>
If a <b>straight pipe with a circular bellmouth inlet (collector) without baffle is mounted flush with the wall</b> then its pressure loss coefficient can be established from the following table. Herein, r is the radius of the bellmouth inlet surface (Idelchik, p. 164 f., Diagram 3-4, paragraph b)
</p>
 
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for outlets, bellmouth flush with wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> r / D_hyd  </th>
  </tr>
  <tr>
    <td></td> <th> 0.01 </th><th> 0.03 </th><th> 0.05 </th><th> 0.08 </th><th> 0.16 </th><th>&ge;0.20</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 0.44 </td><td> 0.31 </td><td> 0.22  </td><td> 0.15  </td><td> 0.06  </td><td>      0.03     </td>
  </tr>
</table>
 
<p>
If a <b>straight pipe with a circular bellmouth inlet (collector) without baffle is mounted at a distance from a wall</b> then its pressure loss coefficient can be established from the following table. Herein, r is the radius of the bellmouth inlet surface (Idelchik, p. 164 f., Diagram 3-4, paragraph a)
</p>
 
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for outlets, bellmouth at a distance of wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> r / D_hyd  </th>
  </tr>
  <tr>
    <td></td> <th> 0.01 </th><th> 0.03 </th><th> 0.05 </th><th> 0.08 </th><th> 0.16 </th><th>&ge;0.20</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 0.87 </td><td> 0.61 </td><td> 0.40  </td><td> 0.20  </td><td> 0.06  </td><td>      0.03     </td>
  </tr>
</table>
 
 
 
<h4><font color=\"#008000\">Inlet Coefficients</font></h4>
 
<p>
If a <b>straight pipe with constant circular cross section is mounted flush with the wall</b>, its vessel inlet pressure loss coefficient will be according to the following table (Idelchik, p. 209 f., Diagram 4-2 with <code>A_port/A_vessel = 0</code> and Idelchik, p. 640, Diagram 11-1, graph a). According to the text, <code>m = 9</code> is appropriate for fully developed turbulent flow.
</p>
 
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for inlets, circular tube flush with wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> m  </th>
  </tr>
  <tr>
    <td></td> <th> 1.0 </th><th> 2.0 </th><th> 3.0 </th><th> 4.0 </th><th> 7.0 </th><th>9.0</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 2.70 </td><td> 1.50 </td><td> 1.25  </td><td> 1.15  </td><td> 1.06  </td><td>      1.04     </td>
  </tr>
</table>
 
 
 
<h4><font color=\"#008000\">References</font></h4>
 
<dl><dt>Idelchik I.E. (1994):</dt>
    <dd><a href=\"http://www.begellhouse.com/books/00c0f05b040d2ec0.html\"><b>Handbook
        of Hydraulic Resistance</b></a>. 3rd edition, Begell House, ISBN
        0-8493-9908-4</dd>
</dl>
</html>"));
    end VesselPortsData;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Vessels;
