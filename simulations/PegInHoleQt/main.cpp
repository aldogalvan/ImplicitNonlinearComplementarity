#include "Interface.h"
#include <QApplication>
//------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
    QApplication a (argc, argv);
    Interface w;
    w.show ();
    w.Start ();

    return a.exec ();
}