//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2021 by Paul Demorest
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Pulsar/DynamicStatSpectrumPlot.h"
#include "Pulsar/DynamicSpectrumPlot.h"
#include "Pulsar/ProfileStatistic.h"

Pulsar::DynamicStatSpectrumPlot::Interface::Interface( 
    DynamicStatSpectrumPlot *instance )
{
  if( instance )
    set_instance( instance );

  add( &DynamicStatSpectrumPlot::get_statistic,
       &DynamicStatSpectrumPlot::set_statistic,
       "stat", "Statistc to evaluate" );

  import( DynamicSpectrumPlot::Interface() );
}
