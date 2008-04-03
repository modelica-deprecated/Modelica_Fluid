within FluidSandbox;
package PressureLosses 
  "Models and functions providing pressure loss correlations " 
  extends Icons.VariantLibrary;
  package BaseClasses 
    extends Icons.BaseClassLibrary;
    partial model PartialWallFriction 
      "Pressure drop in pipe due to wall friction and gravity (for both flow directions)" 
      
      extends Interfaces.PartialComponent;
      
      // Pressure loss function
      replaceable package WallFriction = WallFrictionCorrelations.Laminar 
        extends WallFrictionCorrelations.PartialWallFriction 
        "Characteristic of wall friction"  annotation(choicesAllMatching=true);
      
      // Geometry
      parameter SI.Length length "Length of pipe";
      parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
      parameter SI.Length roughness(min=0) = 2.5e-5 
        "Absolute roughness of pipe (default = smooth steel pipe)" 
          annotation(Dialog(enable=WallFriction.use_roughness));
      
      // Form of pressure loss equation
      parameter Boolean from_dp=true 
        " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
        annotation (Evaluate=true, Dialog(tab="Advanced"));
      
      // Upstream quantities, have to be provided by PartialTwoPortTransport
      Medium.AbsolutePressure p_upstream "Upstream pressure";
      Medium.SpecificEnthalpy h_upstream "Upstream specific mixing enthalpy";
      Medium.MassFraction Xi_upstream[Medium.nXi] "Upstream mass fractions";
      
      // pressure drop quantities, have to be provided by PartialTwoPortTransport
      Medium.MassFlowRate m_flow(start=m_flow_start) 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
      SI.Pressure dp(start=p_a_start - p_b_start) 
        "Pressure drop between port_a and port_b";
      
      //Initialization
      parameter Medium.AbsolutePressure p_a_start 
        "Guess value of pressure at port_a" 
        annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start 
        "Guess value of pressure at port_b" 
        annotation(Dialog(tab = "Initialization"));
      parameter Medium.MassFlowRate m_flow_start 
        "Guesss value of mass flow rate through component" 
        annotation(Dialog(tab = "Initialization"));
      
      // Density and dynamic viscosity for pressure drop correlation
      Medium.Density d_upstream "Upstream density";
      Medium.DynamicViscosity eta_upstream "Upstream dynamics viscosity";
      
    equation 
      d_upstream = Medium.density(Medium.setState_phX(
              p_upstream,
              h_upstream,
              Xi_upstream));
      eta_upstream = Medium.dynamicViscosity(Medium.setState_phX(
              p_upstream,
              h_upstream,
              Xi_upstream));
      
      if from_dp and not WallFriction.dp_is_zero then
        m_flow = WallFriction.massFlowRate_dp(
                dp,
                d_upstream,
                d_upstream,
                eta_upstream,
                eta_upstream,
                length,
                diameter,
                roughness=roughness,
                dp_small=1);
      else
        dp = WallFriction.pressureLoss_m_flow(
                m_flow,
                d_upstream,
                d_upstream,
                eta_upstream,
                eta_upstream,
                length,
                diameter,
                roughness=roughness,
                m_flow_small=0.01);
      end if;
      
      annotation (
        defaultComponentName="pipeFriction",
        Icon(Rectangle(extent=[-100,60; 100,-60], style(
              color=0,
              gradient=2,
              fillColor=8)), Rectangle(extent=[-100,34; 100,-36], style(
              color=69,
              gradient=2,
              fillColor=69))),
        Diagram(
          Rectangle(extent=[-100,64; 100,-64], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=8)),
          Rectangle(extent=[-100,50; 100,-49], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-60,-49; -60,50], style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1)),
          Text(
            extent=[-50,16; 6,-10],
            style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1),
            string="diameter"),
          Line(points=[-100,74; 100,74], style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1)),
          Text(
            extent=[-34,92; 34,74],
            style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1),
            string="length")),
        Coordsys(grid=[1,1], scale=0));
      
    end PartialWallFriction;
    
    partial model PartialWallFrictionWithSmoothing 
      "Pressure drop in pipe due to wall friction and gravity (for both flow directions, with extra smoothing using upstream and downstream conditions)" 
      
      extends Interfaces.PartialComponent;
      
      // Pressure loss function
      replaceable package WallFriction = WallFrictionCorrelations.Laminar 
        extends WallFrictionCorrelations.PartialWallFriction 
        "Characteristic of wall friction"  annotation(choicesAllMatching=true);
      
      // Geometry
      parameter SI.Length length "Length of pipe";
      parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
      parameter SI.Length roughness(min=0) = 2.5e-5 
        "Absolute roughness of pipe (default = smooth steel pipe)" 
          annotation(Dialog(enable=WallFriction.use_roughness));
      
      // Form of pressure loss equation
      parameter Boolean from_dp=true 
        " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
        annotation (Evaluate=true, Dialog(tab="Advanced"));
      
      // pressure drop quantities, have to be provided by PartialTwoPortTransport
      Medium.MassFlowRate m_flow(start=m_flow_start) 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
      SI.Pressure dp(start=p_a_start - p_b_start) 
        "Pressure drop between port_a and port_b";
      
      //Initialization
      parameter Medium.AbsolutePressure p_a_start 
        "Guess value of pressure at port_a" 
        annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start 
        "Guess value of pressure at port_b" 
        annotation(Dialog(tab = "Initialization"));
      parameter Medium.MassFlowRate m_flow_start 
        "Guesss value of mass flow rate through component" 
        annotation(Dialog(tab = "Initialization"));
      
      // Properties for fluid flowing in design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_designDirection 
        "Pressure for flow in design direction";
      Medium.SpecificEnthalpy h_designDirection 
        "Specific mixing enthalpy for flow in design direction";
      Medium.MassFraction Xi_designDirection[Medium.nXi] 
        "Mass fractions for flow in design direction";
      
      // Properties for fluid flowing in non-design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_nonDesignDirection 
        "Pressure for flow in non-design direction";
      Medium.SpecificEnthalpy h_nonDesignDirection 
        "Specific mixing enthalpy for flow in non-design direction";
      Medium.MassFraction Xi_nonDesignDirection[Medium.nXi] 
        "Mass fractions for flow in non-design direction";
      
      // Density and dynamic viscosity for pressure drop correlation
      Medium.Density d_designDirection 
        "Upstream density for design flow direction";
      Medium.DynamicViscosity eta_designDirection 
        "Upstream dynamics viscosity for non-design flow direction";
      Medium.Density d_nonDesignDirection 
        "Upstream density for design flow direction";
      Medium.DynamicViscosity eta_nonDesignDirection 
        "Upstream dynamics viscosity for non-design flow direction";
      
    protected 
      Medium.BaseProperties medium_designDirection(p=p_designDirection, h=h_designDirection, Xi=Xi_designDirection);
      Medium.BaseProperties medium_nonDesignDirection(p=p_nonDesignDirection, h=h_nonDesignDirection, Xi=Xi_nonDesignDirection);
      
    equation 
      d_designDirection = Medium.density(medium_designDirection.state);
      eta_designDirection = if not WallFriction.use_eta then 1.e-10 else Medium.dynamicViscosity(medium_designDirection.state);
      d_nonDesignDirection = Medium.density(medium_nonDesignDirection.state);
      eta_nonDesignDirection = if not WallFriction.use_eta then 1.e-10 else Medium.dynamicViscosity(medium_nonDesignDirection.state);
      
      if from_dp and not WallFriction.dp_is_zero then
        m_flow = WallFriction.massFlowRate_dp(
                dp,
                d_designDirection,
                d_nonDesignDirection,
                eta_designDirection,
                eta_nonDesignDirection,
                length,
                diameter,
                roughness=roughness,
                dp_small=1);
      else
        dp = WallFriction.pressureLoss_m_flow(
                m_flow,
                d_designDirection,
                d_nonDesignDirection,
                eta_designDirection,
                eta_nonDesignDirection,
                length,
                diameter,
                roughness=roughness,
                m_flow_small=0.01);
      end if;
      
      annotation (
        defaultComponentName="pipeFriction",
        Icon(Rectangle(extent=[-100,60; 100,-60], style(
              color=0,
              gradient=2,
              fillColor=8)), Rectangle(extent=[-100,34; 100,-36], style(
              color=69,
              gradient=2,
              fillColor=69))),
        Diagram(
          Rectangle(extent=[-100,64; 100,-64], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=8)),
          Rectangle(extent=[-100,50; 100,-49], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-60,-49; -60,50], style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1)),
          Text(
            extent=[-50,16; 6,-10],
            style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1),
            string="diameter"),
          Line(points=[-100,74; 100,74], style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1)),
          Text(
            extent=[-34,92; 34,74],
            style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1),
            string="length")),
        Coordsys(grid=[1,1], scale=0));
    end PartialWallFrictionWithSmoothing;
  end BaseClasses;
  
  package WallFrictionCorrelations 
    "Different variants for pressure drops due to pipe wall friction (taken one-to-one from Modelica_Fluid)" 
    partial package PartialWallFriction 
      "Partial wall friction characteristic (base package of all wall friction characteristics)" 
      
      annotation (Documentation(info="<html>
 
</html>"));
      
    // Constants to be set in subpackages
      constant Boolean use_eta=true 
        "= true, if eta_a/eta_b are used in function, otherwise value is not used";
      constant Boolean use_roughness=true 
        "= true, if roughness is used in function, otherwise value is not used";
      constant Boolean use_dp_small=true 
        "= true, if dp_small is used in function, otherwise value is not used";
      constant Boolean use_m_flow_small=true 
        "= true, if m_flow_small is used in function, otherwise value is not used";
      constant Boolean dp_is_zero=false 
        "= true, if no wall friction is present, i.e., dp = 0 (function massFlowRate_dp() cannot be used)";
      
    // pressure loss characteristic functions
      replaceable partial function massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        extends Modelica.Icons.Function;
        
        input Modelica.SIunits.Pressure dp 
          "Pressure drop (dp = port_a.p - port_b.p)";
        input Modelica.SIunits.Density d_a "Density at port_a";
        input Modelica.SIunits.Density d_b "Density at port_b";
        input Modelica.SIunits.DynamicViscosity eta_a 
          "Dynamic viscosity at port_a (dummy if use_eta = false)";
        input Modelica.SIunits.DynamicViscosity eta_b 
          "Dynamic viscosity at port_b (dummy if use_eta = false)";
        input Modelica.SIunits.Length length "Length of pipe";
        input Modelica.SIunits.Diameter diameter 
          "Inner (hydraulic) diameter of pipe";
        input Modelica.SIunits.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
        input Modelica.SIunits.AbsolutePressure dp_small=1 
          "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";
        
        output Modelica.SIunits.MassFlowRate m_flow 
          "Mass flow rate from port_a to port_b";
        annotation (Documentation(info="<html>
 
</html>"));
      end massFlowRate_dp;
      
      replaceable partial function pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        extends Modelica.Icons.Function;
        
        input Modelica.SIunits.MassFlowRate m_flow 
          "Mass flow rate from port_a to port_b";
        input Modelica.SIunits.Density d_a "Density at port_a";
        input Modelica.SIunits.Density d_b "Density at port_b";
        input Modelica.SIunits.DynamicViscosity eta_a 
          "Dynamic viscosity at port_a (dummy if use_eta = false)";
        input Modelica.SIunits.DynamicViscosity eta_b 
          "Dynamic viscosity at port_b (dummy if use_eta = false)";
        input Modelica.SIunits.Length length "Length of pipe";
        input Modelica.SIunits.Diameter diameter 
          "Inner (hydraulic) diameter of pipe";
        input Modelica.SIunits.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
        input Modelica.SIunits.MassFlowRate m_flow_small=0.01 
          "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
        output Modelica.SIunits.Pressure dp 
          "Pressure drop (dp = port_a.p - port_b.p)";
        
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
        final use_eta=false,
        final use_roughness=false,
        final use_dp_small=false,
        final use_m_flow_small=false,
        final dp_is_zero=true);
      
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
        final use_eta=true,
        final use_roughness=false,
        final use_dp_small=false,
        final use_m_flow_small=false);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        
        annotation (Documentation(info="<html>
 
</html>"));
      algorithm 
        m_flow := dp*Modelica.Constants.pi*diameter^4*(d_a + d_b)/(128*length
          *(eta_a + eta_b));
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        
        annotation (Documentation(info="<html>
 
</html>"));
      algorithm 
        dp := m_flow*128*length*(eta_a + eta_b)/(Modelica.Constants.pi*
          diameter^4*(d_a + d_b));
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
        final use_eta=false,
        final use_roughness=true,
        final use_dp_small=true,
        final use_m_flow_small=true);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import Modelica.Math;
        import Modelica_Fluid;
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
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
        zeta := (length/diameter)/(2*Math.log10(3.7/(roughness/diameter)))^2;
        k_inv := (pi*diameter*diameter)^2/(8*zeta);
        m_flow := Modelica_Fluid.Utilities.regRoot2(
                  dp,
                  dp_small,
                  d_a*k_inv,
                  d_b*k_inv);
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        import Modelica.Math;
        import Modelica_Fluid;
        
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
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
        zeta := (length/diameter)/(2*Math.log10(3.7/(roughness/diameter)))^2;
        k := 8*zeta/(pi*diameter*diameter)^2;
        dp := Modelica_Fluid.Utilities.regSquare2(
                  m_flow,
                  m_flow_small,
                  k/d_a,
                  k/d_b);
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
        final use_eta=true,
        final use_roughness=true,
        final use_dp_small=false,
        final use_m_flow_small=false);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import Modelica.Math;
        import Modelica_Fluid;
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
        constant Real Re_turbulent=4000 "Start of turbulent regime";
        Real zeta;
        Real k0;
        Real k_inv;
        Real yd0 "Derivative of m_flow=m_flow(dp) at zero";
        Modelica.SIunits.AbsolutePressure dp_turbulent;
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
        zeta := (length/diameter)/(2*Math.log10(3.7/(roughness/diameter)))^2;
        k0 := 128*length/(pi*diameter^4);
        k_inv := (pi*diameter*diameter)^2/(8*zeta);
        yd0 := (d_a + d_b)/(k0*(eta_a + eta_b));
        dp_turbulent := ((eta_a + eta_b)*diameter*pi/8)^2*Re_turbulent^2/(
          k_inv*(d_a + d_b)/2);
        m_flow := Modelica_Fluid.Utilities.regRoot2(
                  dp,
                  dp_turbulent,
                  d_a*k_inv,
                  d_b*k_inv,
                  use_yd0=true,
                  yd0=yd0);
      end massFlowRate_dp;
      
      redeclare function extends pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        import Modelica.Math;
        import Modelica_Fluid;
        
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
        constant Real Re_turbulent=4000 "Start of turbulent regime";
        Real zeta;
        Real k0;
        Real k;
        Real yd0 "Derivative of dp = f(m_flow) at zero";
        Modelica.SIunits.MassFlowRate m_flow_turbulent 
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
        zeta := (length/diameter)/(2*Math.log10(3.7/(roughness/diameter)))^2;
        k0 := 128*length/(pi*diameter^4);
        k := 8*zeta/(pi*diameter*diameter)^2;
        yd0 := k0*(eta_a + eta_b)/(d_a + d_b);
        m_flow_turbulent := (pi/8)*diameter*(eta_a + eta_b)*Re_turbulent;
        dp := Modelica_Fluid.Utilities.regSquare2(
                  m_flow,
                  m_flow_turbulent,
                  k/d_a,
                  k/d_b,
                  use_yd0=true,
                  yd0=yd0);
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
        final use_eta=true,
        final use_roughness=false,
        final use_dp_small=false,
        final use_m_flow_small=false);
      
      redeclare function extends massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import Modelica.Math;
        annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
      protected 
        constant Real pi=Modelica.Constants.pi;
        Real Delta=roughness/diameter "Relative roughness";
        Modelica.SIunits.ReynoldsNumber Re1=(745*Math.exp(if Delta <= 0.0065 then 
                  1 else 0.0065/Delta))^0.97 "Re leaving laminar curve";
        Modelica.SIunits.ReynoldsNumber Re2=4000 "Re entering turbulent curve";
        Modelica.SIunits.DynamicViscosity eta "Upstream viscosity";
        Modelica.SIunits.Density d "Upstream density";
        Modelica.SIunits.ReynoldsNumber Re "Reynolds number";
        Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
        
        function interpolateInRegion2 
          input Real Re_turbulent;
          input Modelica.SIunits.ReynoldsNumber Re1;
          input Modelica.SIunits.ReynoldsNumber Re2;
          input Real Delta;
          input Real lambda2;
          output Modelica.SIunits.ReynoldsNumber Re;
          annotation (smoothOrder=1);
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
        d := if dp >= 0 then d_a else d_b;
        eta := if dp >= 0 then eta_a else eta_b;
        lambda2 := abs(dp)*2*diameter^3*d/(length*eta*eta);
        
        // Determine Re under the assumption of laminar flow
        Re := lambda2/64;
        
        // Modify Re, if turbulent flow
        if Re > Re1 then
          Re := -2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
          if Re < Re2 then
            Re := interpolateInRegion2(
                      Re,
                      Re1,
                      Re2,
                      Delta,
                      lambda2);
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
        constant Real pi=Modelica.Constants.pi;
        Real Delta=roughness/diameter "Relative roughness";
        Modelica.SIunits.ReynoldsNumber Re1=745*Math.exp(if Delta <= 0.0065 then 
                  1 else 0.0065/Delta) "Re leaving laminar curve";
        Modelica.SIunits.ReynoldsNumber Re2=4000 "Re entering turbulent curve";
        Modelica.SIunits.DynamicViscosity eta "Upstream viscosity";
        Modelica.SIunits.Density d "Upstream density";
        Modelica.SIunits.ReynoldsNumber Re "Reynolds number";
        Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
        
        function interpolateInRegion2 
          input Modelica.SIunits.ReynoldsNumber Re;
          input Modelica.SIunits.ReynoldsNumber Re1;
          input Modelica.SIunits.ReynoldsNumber Re2;
          input Real Delta;
          output Real lambda2;
          annotation (smoothOrder=1);
          // point lg(lambda2(Re1)) with derivative at lg(Re1)
        protected 
          Real x1=Math.log10(Re1);
          Real y1=Math.log10(64*Re1);
          Real yd1=1;
          
          // Point lg(lambda2(Re2)) with derivative at lg(Re2)
          Real aux1=(0.5/Math.log(10))*5.74*0.9;
          Real aux2=Delta/3.7 + 5.74/Re2^0.9;
          Real aux3=Math.log10(aux2);
          Real L2=0.25*(Re2/aux3)^2;
          Real aux4=2.51/sqrt(L2) + 0.27*Delta;
          Real aux5=-2*sqrt(L2)*Math.log10(aux4);
          Real x2=Math.log10(Re2);
          Real y2=Math.log10(L2);
          Real yd2=2 + 4*aux1/(aux2*aux3*(Re2)^0.9);
          
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
        d := if m_flow >= 0 then d_a else d_b;
        eta := if m_flow >= 0 then eta_a else eta_b;
        
        // Determine Re, lambda2 and pressure drop
        Re := (4/pi)*abs(m_flow)/(diameter*eta);
        lambda2 := if Re <= Re1 then 64*Re else (if Re >= Re2 then 0.25*(Re/
          Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else interpolateInRegion2(
                  Re,
                  Re1,
                  Re2,
                  Delta));
        dp := length*eta*eta/(2*d*diameter*diameter*diameter)*(if m_flow >= 0 then 
                lambda2 else -lambda2);
      end pressureLoss_m_flow;
    end Detailed;
  end WallFrictionCorrelations;
  
  model WallFriction 
    "Pressure drop in pipe due to wall friction (for both flow directions)" 
    annotation (defaultComponentName="pipeFriction");
    
    extends BaseClasses.PartialWallFrictionWithSmoothing;
    extends FluidInterface.PartialTransportIsenthalpic(
      m_flow(start=m_flow_start),
      port_a(p(start=p_a_start)),
      port_b(p(start=p_b_start)),
      dp(start=p_a_start - p_b_start));
    
  end WallFriction;
  
  model WallFrictionAA 
    "Pressure drop in pipe due to wall friction (for both flow directions, with two PortA's)" 
    annotation (defaultComponentName="pipeFriction");
    
    extends BaseClasses.PartialWallFrictionWithSmoothing;
    extends FluidInterface.PartialTransportIsenthalpicAA(
      m_flow(start=m_flow_start),
      port_a(p(start=p_a_start)),
      port_b(p(start=p_b_start)),
      dp(start=p_a_start - p_b_start));
    
  end WallFrictionAA;
  
  model WallFrictionAB 
    "Pressure drop in pipe due to wall friction (for both flow directions, with a PortA and PortB each)" 
    annotation (defaultComponentName="pipeFriction");
    
    extends BaseClasses.PartialWallFrictionWithSmoothing;
    extends FluidInterface.PartialTransportIsenthalpicAB(
      m_flow(start=m_flow_start),
      port_a(p(start=p_a_start)),
      port_b(p(start=p_b_start)),
      dp(start=p_a_start - p_b_start));
    
  end WallFrictionAB;
end PressureLosses;
