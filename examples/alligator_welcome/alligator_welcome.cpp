#include "tiny_alligator/Alligator.hpp"

int main()
{
    TinyAlligator::play_with_tiny_alligator();
    TinyAlligator::say_bye();
    auto alligator = TinyAlligator::find_tiny_alligator();
    alligator->play_with_tiny_alligator();
    alligator->say_bye();
}