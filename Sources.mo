within Modelica_Fluid;
package Sources
  "Generic sources for fluid connectors to define fixed or prescribed ambient conditions"
  extends Modelica_Fluid.Icons.VariantLibrary;
  model FixedBoundary "Boundary source component"
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean use_p=true "select p or d" 
      annotation (Evaluate = true,
                  Dialog(group = "Boundary pressure or Boundary density"));
    parameter Medium.AbsolutePressure p=Medium.p_default "Boundary pressure" 
      annotation (Dialog(group = "Boundary pressure or Boundary density",
                         enable = use_p));
    parameter Medium.Density d=Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
      "Boundary density" 
      annotation (Dialog(group = "Boundary pressure or Boundary density",
                         enable=not use_p));
    parameter Boolean use_T=true "select T or h" 
      annotation (Evaluate = true,
                  Dialog(group = "Boundary temperature or Boundary specific enthalpy"));
    parameter Medium.Temperature T=Medium.T_default "Boundary temperature" 
      annotation (Dialog(group = "Boundary temperature or Boundary specific enthalpy",
                         enable = use_T));
    parameter Medium.SpecificEnthalpy h=Medium.h_default
      "Boundary specific enthalpy" 
      annotation (Dialog(group="Boundary temperature or Boundary specific enthalpy",
                  enable = not use_T));
    parameter Medium.MassFraction X[Medium.nX](
         quantity=Medium.substanceNames)=Medium.X_default
      "Boundary mass fractions m_i/m" 
      annotation (Dialog(group = "Only for multi-substance flow", enable=Medium.nXi > 0));

    annotation (defaultComponentName = "Boundary_fixed",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={0,127,255}), Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255})}),
      Documentation(info="<html>
<p>
Model <b>FixedBoundary</b> defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure or boundary density.</li>
<li> Boundary temperature or boundary specific enthalpy.</li>
<li> Boundary composition (only for multi-substance flow).</li>
</ul>
<p>
Note, that boundary temperature, density, specific enthalpy
and mass fractions have only an effect if the mass flow
is from the Boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));

  equation
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                          Medium.singleState, use_p, X,
                                          "FixedBoundary");
    if use_p or Medium.singleState then
      medium.p = p;
    else
      medium.d = d;
    end if;
    if use_T then
      medium.T = T;
    else
      medium.h = h;
    end if;

    medium.Xi = X[1:Medium.nXi];
  end FixedBoundary;

  model FixedBoundary_pTX
    "Boundary pressure, temperature and mass fraction source"
    extends Sources.BaseClasses.PartialSource;
    parameter Medium.AbsolutePressure p "Boundary pressure";
    parameter Medium.Temperature T "Boundary temperature";
    parameter Medium.MassFraction X[Medium.nX](
         quantity=Medium.substanceNames) = Medium.X_default
      "Boundary mass fractions m_i/m" 
      annotation (Dialog(group = "Only for multi-substance flow",
                  enable=Medium.nXi > 0));
    annotation (defaultComponentName = "boundary_fixed",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={0,127,255}), Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255})}),
      Documentation(info="<html>
<p>
Defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary temperature.</li>
<li> Boundary composition (only for multi-substance flow).</li>
</ul>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
  equation
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                          Medium.singleState, true, X, "FixedBoundary_pTX");
    medium.p   = p;
    medium.T   = T;
    medium.Xi = X[1:Medium.nXi];
  end FixedBoundary_pTX;

  model FixedBoundary_phX
    "Boundary pressure, specific enthalpy and mass fraction source"
    extends Sources.BaseClasses.PartialSource;
    parameter Medium.AbsolutePressure p "Boundary pressure";
    parameter Medium.SpecificEnthalpy h "Boundary specific enthalpy";
    parameter Medium.MassFraction X[
      Medium.nX](quantity=Medium.substanceNames) = Medium.X_default
      "Boundary mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nXi > 0));
    annotation (defaultComponentName = "boundary_fixed",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={0,127,255}), Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255})}),
      Documentation(info="<html>
<p>
Defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary specific enthalpy.</li>
<li> Boundary composition (only for multi-substance flow).</li>
</ul>
<p>
Note, that boundary specific enthalpy
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
  equation
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                          Medium.singleState, true, X, "FixedBoundary_phX");
    medium.p = p;
    medium.h = h;
    medium.Xi = X[1:Medium.nXi];
  end FixedBoundary_phX;

  annotation (Documentation(info="<html>
<p>
Package <b>Sources</b> contains generic sources for fluid connectors
to define fixed or prescribed ambient conditions.
</p>
</html>"));
  model PrescribedBoundary_pTX
    "Boundary with prescribed pressure, temperature and composition"
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean usePressureInput = false
      "Get the pressure from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useTemperatureInput= false
      "Get the temperature from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.AbsolutePressure p = Medium.p_default
      "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = not usePressureInput));
    parameter Medium.Temperature T = Medium.T_default
      "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = not useTemperatureInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in if              usePressureInput
      "Prescribed boundary pressure" 
      annotation (Placement(transformation(extent={{-140,40},{-100,80}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput T_in if         useTemperatureInput
      "Prescribed boundary temperature" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX] if 
                                                          useCompositionInput
      "Prescribed boundary composition" 
      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}},
            rotation=0)));
  protected
    Modelica.Blocks.Interfaces.RealInput p_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput T_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX]
      "Needed to connect to conditional connector";
    annotation (defaultComponentName = "boundary_prescribed",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={0,127,255}),
          Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255}),
          Line(
            visible=usePressureInput,
            points={{-100,60},{-80,60}},
            color={0,0,255}),
          Line(
            visible=useCompositionInput,
            points={{-100,-60},{-80,-60}},
            color={0,0,255}),
          Text(
            visible=usePressureInput,
            extent={{-146,110},{-62,70}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="p"),
          Text(
            visible=useCompositionInput,
            extent={{-160,-22},{-58,-62}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="X"),
          Text(
            visible=useTemperatureInput,
            extent={{-158,44},{-56,4}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="T")}),
      Documentation(info="<html>
<p>
Defines prescribed values for boundary conditions:
</p>
<ul>
<li> Prescribed boundary pressure.</li>
<li> Prescribed boundary temperature.</li>
<li> Prescribed boundary composition (only for multi-substance flow).</li>
</ul>
<p>If <tt>usePressureInput</tt> is false (default option), the <tt>p</tt> parameter
is used as boundary pressure, and the <tt>p_in</tt> input connector is disabled; if <tt>usePressureInput</tt> is true, then the <tt>p</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics));
  equation
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedBoundary_pTX");
    connect(p_in, p_in_internal);
    connect(T_in, T_in_internal);
    connect(X_in, X_in_internal);
    if not usePressureInput then
      p_in_internal = p;
    end if;
    if not useTemperatureInput then
      T_in_internal = T;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    medium.p = p_in_internal;
    medium.T = T_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedBoundary_pTX;

  model PrescribedBoundary_phX
    "Boundary with prescribed pressure, specific enthalpy and composition"
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean usePressureInput = false
      "Get the pressure from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useEnthalpyInput= false
      "Get the specific enthalpy from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.AbsolutePressure p = Medium.p_default
      "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = not usePressureInput));
    parameter Medium.SpecificEnthalpy h = Medium.h_default
      "Fixed value of specific enthalpy" 
      annotation (Evaluate = true,
                  Dialog(enable = not useEnthalpyInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in if              usePressureInput
      "Prescribed boundary pressure" 
      annotation (Placement(transformation(extent={{-140,40},{-100,80}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput h_in if              useEnthalpyInput
      "Prescribed boundary specific enthalpy" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX] if 
                                                          useCompositionInput
      "Prescribed boundary composition" 
      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}},
            rotation=0)));
    annotation (defaultComponentName = "boundary_prescribed",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={0,127,255}),
          Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255}),
          Line(
            visible=usePressureInput,
            points={{-100,60},{-66,60}},
            color={0,0,255}),
          Line(
            visible=useCompositionInput,
            points={{-100,-60},{-66,-60}},
            color={0,0,255}),
          Text(
            visible=usePressureInput,
            extent={{-148,120},{-70,80}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="p"),
          Text(
            visible=useEnthalpyInput,
            extent={{-100,20},{2,-20}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="h"),
          Text(
            visible=useCompositionInput,
            extent={{-140,-86},{-38,-126}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="X")}),
      Documentation(info="<html>
<p>
Defines prescribed values for boundary conditions:
</p>
<ul>
<li> Prescribed boundary pressure.</li>
<li> Prescribed boundary temperature.</li>
<li> Prescribed boundary composition (only for multi-substance flow).</li>
</ul>
<p>If <tt>usePressureInput</tt> is false (default option), the <tt>p</tt> parameter
is used as boundary pressure, and the <tt>p_in</tt> input connector is disabled; if <tt>usePressureInput</tt> is true, then the <tt>p</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the specific enthalpy and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
  protected
    Modelica.Blocks.Interfaces.RealInput p_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput h_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX]
      "Needed to connect to conditional connector";
  equation
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedBoundary_phX");
    connect(p_in, p_in_internal);
    connect(h_in, h_in_internal);
    connect(X_in, X_in_internal);
    if not usePressureInput then
      p_in_internal = p;
    end if;
    if not useEnthalpyInput then
      h_in_internal = h;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    medium.p = p_in_internal;
    medium.h = h_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedBoundary_phX;

  model PrescribedMassFlowRate_TX
    "Ideal flow source that produces a prescribed mass flow with prescribed temperature and mass fraction"
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean useFlowRateInput = false
      "Get the mass flow rate from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useTemperatureInput= false
      "Get the temperature from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.MassFlowRate m_flow = 0
      "Fixed mass flow rate going out of the fluid port" 
      annotation (Evaluate = true,
                  Dialog(enable = not useFlowRateInput));
    parameter Medium.Temperature T = Medium.T_default
      "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = not useTemperatureInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in if     useFlowRateInput
      "Prescribed mass flow rate" 
      annotation (Placement(transformation(extent={{-113,40},{-73,80}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput T_in if         useTemperatureInput
      "Prescribed fluid temperature" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX] if 
                                                          useCompositionInput
      "Prescribed fluid composition" 
      annotation (Placement(transformation(extent={{-112,-81},{-72,-41}},
            rotation=0)));
  protected
    Modelica.Blocks.Interfaces.RealInput m_flow_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput T_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX]
      "Needed to connect to conditional connector";
    annotation (defaultComponentName = "massFlowRate",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{20,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{38,40},{100,-40}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Ellipse(
            extent={{-100,80},{60,-80}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-60,70},{60,0},{-60,-68},{-60,70}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-54,32},{16,-30}},
            lineColor={255,0,0},
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid,
            textString="m"),
          Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255}),
          Ellipse(
            extent={{-26,30},{-18,22}},
            lineColor={255,0,0},
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            visible=useFlowRateInput,
            extent={{-194,112},{-54,80}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="m_flow"),
          Text(
            visible=useTemperatureInput,
            extent={{-100,14},{-60,-20}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="T"),
          Text(
            visible=useCompositionInput,
            extent={{-144,-90},{-24,-118}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="X")}),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
      Documentation(info="<html>
<p>
Models an ideal flow source, with prescribed values of flow rate, temperature and composition:
</p>
<ul>
<li> Prescribed mass flow rate.</li>
<li> Prescribed temperature.</li>
<li> Prescribed composition (only for multi-substance flow) .</li>
</ul>
<p>If <tt>useFlowRateInput</tt> is false (default option), the <tt>m_flow</tt> parameter
is used as boundary pressure, and the <tt>m_flow_in</tt> input connector is disabled; if <tt>useFlowRateInput</tt> is true, then the <tt>m_flow</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary flow rate, do not have an effect.
</p>
</html>"));
  equation
    Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedMassFlowRate_TX");
    connect(m_flow_in, m_flow_in_internal);
    connect(T_in, T_in_internal);
    connect(X_in, X_in_internal);
    if not useFlowRateInput then
      m_flow_in_internal = m_flow;
    end if;
    if not useTemperatureInput then
      T_in_internal = T;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    port.m_flow = -m_flow_in_internal;
    medium.T = T_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedMassFlowRate_TX;

  model PrescribedMassFlowRate_hX
    "Ideal flow source that produces a prescribed mass flow with prescribed specific enthalpy and mass fraction"
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean useFlowRateInput = false
      "Get the mass flow rate from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useEnthalpyInput= false
      "Get the specific enthalpy from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.MassFlowRate m_flow = 0
      "Fixed mass flow rate going out of the fluid port" 
      annotation (Evaluate = true,
                  Dialog(enable = not useFlowRateInput));
    parameter Medium.SpecificEnthalpy h = Medium.h_default
      "Fixed value of specific enthalpy" 
      annotation (Evaluate = true,
                  Dialog(enable = not useEnthalpyInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in if     useFlowRateInput
      "Prescribed mass flow rate" 
      annotation (Placement(transformation(extent={{-113,40},{-73,80}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput h_in if              useEnthalpyInput
      "Prescribed fluid specific enthalpy" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX] if 
                                                          useCompositionInput
      "Prescribed fluid composition" 
      annotation (Placement(transformation(extent={{-113,-80},{-73,-40}},
            rotation=0)));
  protected
    Modelica.Blocks.Interfaces.RealInput m_flow_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput h_in_internal
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX]
      "Needed to connect to conditional connector";
    annotation (defaultComponentName = "massFlowRate",
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{20,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{38,40},{100,-40}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Ellipse(
            extent={{-100,80},{60,-80}},
            lineColor={0,0,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-60,70},{60,0},{-60,-68},{-60,70}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-54,32},{16,-30}},
            lineColor={255,0,0},
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid,
            textString="m"),
          Ellipse(
            extent={{-26,30},{-18,22}},
            lineColor={255,0,0},
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            visible=useFlowRateInput,
            extent={{-194,115},{-54,83}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="m_flow"),
          Text(
            visible=useEnthalpyInput,
            extent={{-100,15},{-60,-19}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="h"),
          Text(
            visible=useCompositionInput,
            extent={{-145,-85},{-25,-113}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="X"),
          Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255})}),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
      Documentation(info="<html>
<p>
Models an ideal flow source, with prescribed values of flow rate, temperature and composition:
</p>
<ul>
<li> Prescribed mass flow rate.</li>
<li> Prescribed specific enthalpy.</li>
<li> Prescribed composition (only for multi-substance flow) .</li>
</ul>
<p>If <tt>useFlowRateInput</tt> is false (default option), the <tt>m_flow</tt> parameter
is used as boundary pressure, and the <tt>m_flow_in</tt> input connector is disabled; if <tt>useFlowRateInput</tt> is true, then the <tt>m_flow</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary flow rate, do not have an effect.
</p>
</html>"));
  equation
    Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedMassFlowRate_hX");
    connect(m_flow_in, m_flow_in_internal);
    connect(h_in, h_in_internal);
    connect(X_in, X_in_internal);
    if not useFlowRateInput then
      m_flow_in_internal = m_flow;
    end if;
    if not useEnthalpyInput then
      h_in_internal = h;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    port.m_flow = -m_flow_in_internal;
    medium.h = h_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedMassFlowRate_hX;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;
  partial model PartialSource
      "Partial component source with one fluid connector"
      import Modelica.Constants;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium
        "Medium model within the source" 
       annotation (choicesAllMatching=true);

    Medium.BaseProperties medium "Medium in the source";
    Modelica_Fluid.Interfaces.FluidPort_b port(
                                redeclare package Medium = Medium,
                     m_flow(max=if flowDirection==Types.PortFlowDirection.Leaving then 0 else 
                                     +Constants.inf,
                            min=if flowDirection==Types.PortFlowDirection.Entering then 0 else 
                                     -Constants.inf)) 
      annotation (Placement(transformation(extent={{90,-10},{110,10}}, rotation=0)));
    annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_flow) and the mass flow rates of the
     substances (= port.mX_flow) depend on the direction of the mass flow rate.</li>
</ul>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
    protected
    parameter Types.PortFlowDirection flowDirection=
                     Types.PortFlowDirection.Bidirectional
        "Allowed flow direction"             annotation(Evaluate=true, Dialog(tab="Advanced"));
  equation
    port.p = medium.p;
    port.h_outflow  = medium.h;
    port.Xi_outflow = medium.Xi;
  end PartialSource;
  end BaseClasses;
end Sources;
