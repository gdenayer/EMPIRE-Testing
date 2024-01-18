#include "cppunit/extensions/TestFactoryRegistry.h"
#include "cppunit/ui/text/TestRunner.h"
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

class Intt {
private:
    int number;
public:
    Intt(int _number) :
            number(_number) {
    }
    int getNumber() {
        return number;
    }
};

class InttTest;
CPPUNIT_TEST_SUITE_REGISTRATION( InttTest );
class InttTest: public CppUnit::TestFixture {
private:
    Intt *one;
    Intt *two;
public:
    void setUp() {
        one = new Intt(1);
        two = new Intt(2);
    }
    void tearDown() {
        delete one;
        delete two;
    }
    void testOne() {
        CPPUNIT_ASSERT(one->getNumber() == 1);
    }

    void testTwo() {
        CPPUNIT_ASSERT(two->getNumber() == 2);
    }
CPPUNIT_TEST_SUITE( InttTest );
        CPPUNIT_TEST( testOne);
        CPPUNIT_TEST( testTwo);
    CPPUNIT_TEST_SUITE_END();
};


int main(int argc, char **argv) {
    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest(registry.makeTest());
    runner.run();
    return 0;
}
