within FluidSandbox;
package Sources "Source models" 
  extends Icons.VariantLibrary;
  
  model PrescribedBoundary_pTX_A 
    "Boundary with prescribed pressure, temperature and composition" 
    
    extends BaseClasses.PartialPrescribedBoundary_pTX;
    extends FluidInterface.PartialSource_A;
    
  end PrescribedBoundary_pTX_A;
  
  model PrescribedBoundary_pTX_B 
    "Boundary with prescribed pressure, temperature and composition" 
    
    extends BaseClasses.PartialPrescribedBoundary_pTX;
    extends FluidInterface.PartialSource_B;
    
  end PrescribedBoundary_pTX_B;
  
  model PrescribedMassFlowRate_TX_A 
    "Ideal flow source that produces a prescribed mass flow with prescribed temperature and mass fraction" 
    
    extends BaseClasses.PartialPrescribedMassFlowRate_TX;
    extends FluidInterface.PartialSource_A;
    
  end PrescribedMassFlowRate_TX_A;
  
  model PrescribedMassFlowRate_TX_B 
    "Ideal flow source that produces a prescribed mass flow with prescribed temperature and mass fraction" 
    
    extends BaseClasses.PartialPrescribedMassFlowRate_TX;
    extends FluidInterface.PartialSource_B;
    
  end PrescribedMassFlowRate_TX_B;
  
  model ControlledPump "Simple model of pump with prescribed mass flow rate" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsentropic;
    
    Modelica.Blocks.Interfaces.RealInput massFlowRate(redeclare type SignalType
        =            SI.MassFlowRate) 
      annotation (extent=[-10,80; 10,100],  rotation=270);
  equation 
    port_a.m_flow = massFlowRate;
    annotation (Icon(
        Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64], style(
              pattern=0, fillColor=74)),
        Ellipse(extent=[-80,80; 80,-80], style(gradient=3)),
        Polygon(points=[-38,40; -38,-40; 54,0; -38,40], style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-120; 100,-140], string="%name"),
        Text(extent=[20,100; 48,80], string="m_flow"),
        Line(points=[-100,0; -80,0], style(color=3, rgbcolor={0,0,255})),
        Line(points=[80,0; 100,0], style(color=3, rgbcolor={0,0,255}))), Diagram);
  end ControlledPump;
  
  model ControlledPumpAA 
    "Simple model of pump with prescribed mass flow rate (two PortA's, not supported for all interfaces)" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsentropicAA;
    
    Modelica.Blocks.Interfaces.RealInput m_flow(redeclare type SignalType = 
                     SI.MassFlowRate) 
      annotation (extent=[-10,80; 10,100],  rotation=270);
  equation 
    port_a.m_flow = m_flow;
    annotation (Icon(
        Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64], style(
              pattern=0, fillColor=74)),
        Ellipse(extent=[-80,80; 80,-80], style(gradient=3)),
        Polygon(points=[-38,40; -38,-40; 54,0; -38,40], style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-120; 100,-140], string="%name"),
        Text(extent=[20,100; 48,80], string="m_flow"),
        Line(points=[-100,0; -80,0], style(color=3, rgbcolor={0,0,255})),
        Line(points=[80,0; 100,0], style(color=3, rgbcolor={0,0,255}))), Diagram);
  end ControlledPumpAA;

  model ControlledPumpAB 
    "Simple model of pump with prescribed mass flow rate (PortA and PortB, not supported for all interfaces)" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsentropicAB;
    
    Modelica.Blocks.Interfaces.RealInput m_flow(redeclare type SignalType = 
                     SI.MassFlowRate) 
      annotation (extent=[-10,80; 10,100],  rotation=270);
  equation 
    port_a.m_flow = m_flow;
    annotation (Icon(
        Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64], style(
              pattern=0, fillColor=74)),
        Ellipse(extent=[-80,80; 80,-80], style(gradient=3)),
        Polygon(points=[-38,40; -38,-40; 54,0; -38,40], style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-120; 100,-140], string="%name"),
        Text(extent=[20,100; 48,80], string="m_flow"),
        Line(points=[-100,0; -80,0], style(color=3, rgbcolor={0,0,255})),
        Line(points=[80,0; 100,0], style(color=3, rgbcolor={0,0,255}))), Diagram);
  end ControlledPumpAB;

  package BaseClasses 
    extends Icons.BaseClassLibrary;
    partial model PartialPrescribedBoundary_pTX 
      "Partial boundary with prescribed pressure, temperature and composition" 
      
      extends Interfaces.PartialComponent;
      
      // Medium instance  
      Medium.BaseProperties medium "Medium in the source";
      
      // Flags to use conditional input connectors
      parameter Boolean usePressureInput=false 
        "Get the pressure from the input connector";
      parameter Boolean useTemperatureInput=false 
        "Get the temperature from the input connector";
      parameter Boolean useCompositionInput=false 
        "Get the composition from the input connector";
      
      // Source parameters in user interface
      parameter Medium.AbsolutePressure p=Medium.reference_p 
        "Fixed value of pressure" annotation (Evaluate=true,
               Dialog(enable = not usePressureInput));
      parameter Medium.Temperature T=Medium.reference_T 
        "Fixed value of temperature" annotation (Evaluate=true,
               Dialog(enable = not useTemperatureInput));
      parameter Medium.MassFraction X[Medium.nX]=Medium.X_default 
        "Fixed value of composition" annotation (Evaluate=true,
               Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
      
      // Source property input connectors
      Modelica.Blocks.Interfaces.RealInput p_in(redeclare type SignalType = 
            Medium.AbsolutePressure) if                    usePressureInput 
        "Prescribed boundary pressure" annotation (extent=[-140,40; -100,80]);
      Modelica.Blocks.Interfaces.RealInput T_in(redeclare type SignalType = 
            Medium.Temperature) if                    useTemperatureInput 
        "Prescribed boundary temperature" annotation (extent=[-140,-20; -100,20]);
      Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](redeclare type SignalType = 
            Medium.MassFraction) if                    useCompositionInput 
        "Prescribed boundary composition" annotation (extent=[-140,-80; -100,-40]);
      
      // Interface to non-generic (i.e. approach-specific) implementations
    protected 
      Medium.MassFlowRate port_m_flow "Port mass flow rate";
      
      // Internal signal routing
      Modelica.Blocks.Interfaces.RealInput p_in_internal(redeclare type 
          SignalType = Medium.AbsolutePressure) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput T_in_internal(redeclare type 
          SignalType = Medium.Temperature) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](redeclare 
          type SignalType = Medium.MassFraction) 
        "Needed to connect to conditional connector";
      annotation (
        defaultComponentName="boundary_prescribed", 
        Coordsys(
          extent=[-100,-100; 100,100], 
          grid=[2,2], 
          component=[20,20]), 
        Icon(
          Ellipse(extent=[-100,100; 100,-100], style(
              color=69, 
              gradient=3, 
              fillColor=69)), 
          Line(points=[-100,60; -80,60], style(color=3, rgbcolor={0,0,255})), 
          Line(points=[-100,-60; -80,-60], style(color=3, rgbcolor={0,0,255})), 
            
          Text(
            extent=[-146,110; -62,70], 
            style(
              color=0, 
              rgbcolor={0,0,0}, 
              fillColor=7, 
              rgbfillColor={255,255,255}), 
            string="p"), 
          Text(
            extent=[-160,-22; -58,-62], 
            style(
              color=0, 
              rgbcolor={0,0,0}, 
              fillColor=7, 
              rgbfillColor={255,255,255}), 
            string="X"), 
          Text(
            extent=[-158,44; -56,4], 
            style(
              color=0, 
              rgbcolor={0,0,0}, 
              fillColor=7, 
              rgbfillColor={255,255,255}), 
            string="T")), 
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
        Diagram);
      
    equation 
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
    end PartialPrescribedBoundary_pTX;
    
    partial model PartialPrescribedMassFlowRate_TX 
      "Partial ideal flow source that produces a prescribed mass flow with prescribed temperature and mass fraction" 
      
      extends Interfaces.PartialComponent;
      
      // Medium instance  
      Medium.BaseProperties medium "Medium in the source";
      
      // Flags to use conditional input connectors
      parameter Boolean useFlowRateInput=false 
        "Get the mass flow rate from the input connector";
      parameter Boolean useTemperatureInput=false 
        "Get the temperature from the input connector";
      parameter Boolean useCompositionInput=false 
        "Get the composition from the input connector";
      
      // Source parameters in user interface
      parameter Medium.MassFlowRate m_flow=0 
        "Fixed mass flow rate going out of the fluid port" 
      annotation (Evaluate = true,
                  Dialog(enable = not useFlowRateInput));
      parameter Medium.Temperature T=Medium.reference_T 
        "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = not useTemperatureInput));
      parameter Medium.MassFraction X[Medium.nX]=Medium.X_default 
        "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
      
      // Source property input connectors
      Modelica.Blocks.Interfaces.RealInput m_flow_in(redeclare type SignalType 
          = Medium.MassFlowRate) if                       useFlowRateInput 
        "Prescribed mass flow rate" 
      annotation (extent=[-113,40; -73,80]);
      Modelica.Blocks.Interfaces.RealInput T_in(redeclare type SignalType = 
            Medium.Temperature) if                       useTemperatureInput 
        "Prescribed fluid temperature" 
      annotation (extent=[-140,-20; -100,20]);
      Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](redeclare type 
          SignalType = Medium.MassFraction) if            useCompositionInput 
        "Prescribed fluid composition" 
      annotation (extent=[-112,-81; -72,-41]);
      
      // Interface to non-generic (i.e. approach-specific) implementations
    protected 
      Medium.MassFlowRate port_m_flow "Port mass flow rate";
      
      // Internal signal routing
      Modelica.Blocks.Interfaces.RealInput m_flow_in_internal(redeclare type 
          SignalType = Medium.MassFlowRate) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput T_in_internal(redeclare type 
          SignalType = Medium.Temperature) 
        "Needed to connect to conditional connector";
      Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](redeclare 
          type SignalType = Medium.MassFraction) 
        "Needed to connect to conditional connector";
      annotation (
        defaultComponentName="massFlowRate", 
        Coordsys(
          extent=[-100,-100; 100,100], 
          grid=[2,2], 
          component=[20,20], 
          scale=0), 
        Icon(
          Rectangle(extent=[20,60; 100,-60], style(
              color=0, 
              gradient=2, 
              fillColor=8)), 
          Rectangle(extent=[38,40; 100,-40], style(
              color=69, 
              gradient=2, 
              fillColor=69)), 
          Ellipse(extent=[-100,80; 60,-80], style(fillColor=7)), 
          Polygon(points=[-60,70; 60,0; -60,-68; -60,70], style(color=73, 
                fillColor=73)), 
          Text(
            extent=[-54,32; 16,-30], 
            style(color=41, fillColor=41), 
            string="m"), 
          Ellipse(extent=[-26,30; -18,22], style(color=1, fillColor=1)), 
          Text(
            extent=[-194,112; -54,80], 
            style(
              color=0, 
              rgbcolor={0,0,0}, 
              fillColor=7, 
              rgbfillColor={255,255,255}), 
            string="m_flow"), 
          Text(
            extent=[-100,14; -60,-20], 
            style(
              color=0, 
              rgbcolor={0,0,0}, 
              fillColor=7, 
              rgbfillColor={255,255,255}), 
            string="T"), 
          Text(
            extent=[-144,-90; -24,-118], 
            style(
              color=0, 
              rgbcolor={0,0,0}, 
              fillColor=7, 
              rgbfillColor={255,255,255}), 
            string="X")), 
        Window(
          x=0.45, 
          y=0.01, 
          width=0.44, 
          height=0.65), 
        Diagram, 
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
      port_m_flow = -m_flow_in_internal;
      medium.T = T_in_internal;
      medium.Xi = X_in_internal[1:Medium.nXi];
    end PartialPrescribedMassFlowRate_TX;
  end BaseClasses;
end Sources;
