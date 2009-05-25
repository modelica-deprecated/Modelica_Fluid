within Modelica_Fluid.Test.TestComponents.Sources;
model TestSources01

  Modelica_Fluid.Sources.Boundary_pT sink(nPorts=1, redeclare package Medium = 
        ModelicaNew.Media.FluidTypes.IdealGases.NasaAirN2O2) 
    annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
  Modelica_Fluid.Sources.MassFlowSource_T source(
    nPorts=1,
    redeclare package Medium = 
        ModelicaNew.Media.FluidTypes.IdealGases.NasaAirN2O2,
    m_flow=1) annotation (Placement(transformation(extent={{-20,10},{0,30}})));
equation
  connect(source.ports[1], sink.ports[1]) annotation (Line(
      points={{0,20},{12,20},{12,-20},{0,-20}},
      color={0,127,255},
      smooth=Smooth.None));
end TestSources01;
