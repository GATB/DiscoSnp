

// We include what we need for the test
#include <Kissreads2.hpp>
//#include <gatb/gatb_core.hpp>
using namespace std;

/********************************************************************************/
// Once our tool class is defined, we can run it in the main function of the program.
int main (int argc, char* argv[])
{
    // We use a try/catch block since GATB functions may throw exceptions
    try
    {
        // We run our tool with the provided command line arguments.
        // This will call the Kissreads2::execute method we have defined.
        Kissreads2().run (argc, argv);
        // You can try to launch our tool with different command line arguments.
        // For instance, you can try different number of threads:
        // ./ToyTool -max 200000 -nb-cores 1
        // ./ToyTool -max 200000 -nb-cores 2
        // ./ToyTool -max 200000 -nb-cores 4
        // ./ToyTool -max 200000 -nb-cores 8
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}


