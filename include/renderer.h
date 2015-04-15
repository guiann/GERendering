#ifndef RENDERER_H
#define RENDERER_H


class renderer
{
    public:
        renderer();
        virtual ~renderer();
        unsigned int Gettest() { return m_test; }
        void Settest(unsigned int val) { m_test = val; }
    protected:
    private:
        unsigned int m_test;
};

#endif // RENDERER_H
