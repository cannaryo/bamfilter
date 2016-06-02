#include "sequence-splitter.h"
#include "bamfilter-main.h"

int main(int argc, char *argv[])
{
  fmu_tools::BamFilterMain app;
  return app.Run(argc, argv);
}
