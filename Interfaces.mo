within Modelica_Fluid;
package Interfaces
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow"

  annotation (Documentation(info="<html>
 
</html>", revisions="<html>
<ul>
<li><i>June 9th, 2008</i>
       by Michael Sielemann: Introduced stream keyword after decision at 57th Design Meeting (Lund).</li>
<li><i>May 30, 2007</i>
       by Christoph Richter: moved everything back to its original position in Modelica_Fluid.</li>
<li><i>Apr. 20, 2007</i>
       by Christoph Richter: moved parts of the original package from Modelica_Fluid
       to the development branch of Modelica 2.2.2.</li>
<li><i>Nov. 2, 2005</i>
       by Francesco Casella: restructured after 45th Design Meeting.</li>
<li><i>Nov. 20-21, 2002</i>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.
<li><i>Nov. 11, 2002</i>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 6, 2002</i>
       by Hilding Elmqvist: first version.</li>
<li><i>Aug. 11, 2002</i>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</li>
</ul>
</html>"));

  extends Modelica.Icons.Library;

  connector FluidPort
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"

    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
      "Medium model" annotation (choicesAllMatching=true);

    flow Medium.MassFlowRate m_flow
      "Mass flow rate from the connection point into the component";
    Medium.AbsolutePressure p "Pressure in the connection point";
    stream Medium.SpecificEnthalpy h_outflow
      "Specific enthalpy close to the connection point if m_flow < 0";
    stream Medium.MassFraction Xi_outflow[Medium.nXi]
      "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
    stream Medium.ExtraProperty C_outflow[Medium.nC]
      "Properties c_i/m close to the connection point if m_flow < 0";
  end FluidPort;

  connector FluidPort_a "Generic fluid connector at design inlet"
    extends FluidPort;
    annotation (defaultComponentName="port_a",
                Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid), Text(extent={{-150,110},{150,50}},
              textString="%name")}),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,127,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPort_a;

  connector FluidPort_b "Generic fluid connector at design outlet"
    extends FluidPort;
    annotation (defaultComponentName="port_b",
                Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(extent={{-150,110},{150,50}}, textString="%name")}),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,127,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPort_b;

  connector FluidStatePort_a
    "Fluid connector at design inlet with potential pressure state"
    extends FluidPort;
    annotation (defaultComponentName="port_a",
                Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={0,0,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid), Text(extent={{-150,110},{150,50}},
              textString="%name")}),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,127,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-18,20},{22,-20}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end FluidStatePort_a;

  connector FluidStatePort_b
    "Fluid connector at design outlet with potential pressure state"
   extends FluidPort;
    annotation (defaultComponentName="port_b",
                Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={0,0,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(extent={{-150,110},{150,50}}, textString="%name")}),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,127,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-18,20},{22,-20}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end FluidStatePort_b;

  connector FluidPorts_a
    "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
    extends FluidPort;
    annotation (defaultComponentName="ports_a",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,90},{25,40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,25},{25,-25}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,-40},{25,-90}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,180},{50,80}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,-80},{50,-180}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPorts_a;

  connector FluidPorts_b
    "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
    extends FluidPort;
    annotation (defaultComponentName="ports_b",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,90},{25,40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,25},{25,-25}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,-40},{25,-90}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,-50},{15,-80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,15},{15,-15}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,50},{15,80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,180},{50,80}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,-80},{50,-180}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,100},{30,160}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,-100},{30,-160}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPorts_b;

  connector FluidStatePorts_a
    "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
    extends FluidPort;
    annotation (defaultComponentName="ports_a",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,90},{25,40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,25},{25,-25}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,-40},{25,-90}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,180},{50,80}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,-80},{50,-180}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-20,150},{20,110}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-20,20},{20,-20}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-19,-111},{21,-151}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end FluidStatePorts_a;

  connector FluidStatePorts_b
    "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
    extends FluidPort;
    annotation (defaultComponentName="ports_b",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,90},{25,40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,25},{25,-25}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,-40},{25,-90}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,-50},{15,-80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,15},{15,-15}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,50},{15,80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,180},{50,80}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,-80},{50,-180}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-45,44},{45,-44}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-43,88},{43,172}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-43,-84},{45,-175}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-20,151},{20,111}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-20,21},{20,-19}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-19,-110},{21,-150}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end FluidStatePorts_b;

  connector HeatPorts_a
    "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)"
    extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    annotation (defaultComponentName="heatPorts_a",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString=               "%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-21,89},{22,44}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-21,22},{22,-23}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-21.5,-43},{21.5,-88}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-44,176},{44,86}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-43,46},{45,-44}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-44,-86},{44,-176}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid)}));
  end HeatPorts_a;

  connector HeatPorts_b
    "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)"
    extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    annotation (defaultComponentName="heatPorts_b",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-21,88},{22,43}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-21,22},{22,-23}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-21,-43},{22,-88}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-43,175},{45,85}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-44,46},{44,-44}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-44,-86},{44,-176}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end HeatPorts_b;
end Interfaces;
